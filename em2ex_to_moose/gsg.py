"""
gsg.py - reader for Petrel's binary "GSG" simulation grid (the file a Petrel-exported Eclipse
deck loads with the PETGRID keyword when GRIDFILE does not ask for an EGRID).

The layout below was worked out from a Petrel 2017.1 export and checked cell-by-cell against
the Eclipse INIT file of the same run (TOPS/DEPTH/DZ to < 2 mm on every active cell, and the
property arrays bit-identical). Treat other Petrel versions with care: `read_gsg` returns
self-check numbers you should look at (the converter prints them).

File = sequence of tagged sections + an INDEX at the end (little-endian):
  string  := int32 n, n bytes
  INDEX   := int32 0, int32 n, then n x (string tag, int32 0, int64 offset-of-section-body)
  AXES    := 3 unit strings, float32 x0, y0, ..., (x/y map origin of local coordinates)
  GRID    := int32 0, string name, int32 -1, int32 nx, ny, nz, ...
  AREAL   := int32 -1, 0, ncol, ?  then ncol x [col, 4, p0, p1, p2, p3]
             pillars of column (i,r): p0=(i,r+1) p1=(i,r) p2=(i+1,r) p3=(i+1,r+1), id = r*(nx+1)+i
  PILLARS := int32 0, npil, nlev, 0,0,-1,0, 9 ints, then npil records (trailing empties omitted):
               -1                                   empty pillar
               2, f32[6] line(x1,y1,z1,x2,y2,z2), int32 nseq, [int32 slot->seq map[4] if nseq>1],
               f32[nseq*nlev] elevations (FLT_MAX = undefined)
             the column at AREAL corner q uses seq map[3-q]
  PROP    := int32 0, int32 idx, string name, 'c   ', int32 isfloat, 0, count, data(count)
GSG rows run opposite to Eclipse J (eclipse j = ny-1-r); z is elevation in metres; x,y are
metres relative to the AXES origin.
"""
from __future__ import annotations

import re
import struct

import numpy as np

FLT_MAX_CUT = 1e37


def _index(d):
    """Return {tag: [offsets]} from the INDEX table at the end of the file."""
    pos = d.rfind(b"INDEX")
    if pos < 0:
        raise ValueError("not a GSG file (no INDEX)")
    p = pos + 5 + 8  # tag, int32 0, int32 number of entries
    out = {}
    while p + 4 <= len(d):
        (n,) = struct.unpack_from("<i", d, p)
        if n <= 0 or n > 32 or p + 4 + n + 12 > len(d):
            break
        tag = d[p + 4:p + 4 + n].decode("latin1")
        (off,) = struct.unpack_from("<q", d, p + 4 + n + 4)
        out.setdefault(tag, []).append(off)
        p += 4 + n + 12
    return out


def read_gsg(path, read_props=True):
    d = open(path, "rb").read()
    if not d.startswith(b"GSG"):
        raise ValueError("%s is not a GSG file" % path)
    idx = _index(d)

    # AXES: origin of local x,y
    a = idx["AXES"][0]
    # three 4-char unit strings preceded by an int32 0, then int32 4?, then float32 x0, y0
    ax = np.frombuffer(d[a:a + 64], "<f4")
    units = d[a:a + 12].decode("latin1")
    x0, y0 = float(ax[5]), float(ax[6])

    g = idx["GRID"][0]
    nl = struct.unpack_from("<i", d, g + 4)[0]
    name = d[g + 8:g + 8 + nl].decode("latin1")
    q = g + 8 + nl
    nx, ny, nz = [int(v) for v in np.frombuffer(d[q + 4:q + 16], "<i4")]

    ar = idx["AREAL"][0]
    hdr = np.frombuffer(d[ar:ar + 16], "<i4")
    ncol = int(hdr[2])
    if ncol != nx * ny:
        raise ValueError("AREAL column count %d != nx*ny %d" % (ncol, nx * ny))
    rec = np.frombuffer(d[ar + 16:ar + 16 + 24 * ncol], "<i4").reshape(ncol, 6)
    if not (rec[:, 0] == np.arange(ncol)).all() or not (rec[:, 1] == 4).all():
        raise ValueError("unexpected AREAL layout")
    colp = rec[:, 2:6].copy()

    pl = idx["PILLARS"][0]
    end = min([o for t, offs in idx.items() for o in offs if o > pl] + [len(d)])
    I = np.frombuffer(d[pl:pl + (end - pl) // 4 * 4], "<i4")
    F = I.view("<f4")
    npil, nlev = int(I[1]), int(I[2])
    if npil != (nx + 1) * (ny + 1) or nlev != nz + 1:
        raise ValueError("PILLARS header %d/%d does not match grid %dx%dx%d" % (npil, nlev, nx, ny, nz))
    k = 16
    lines = np.full((npil, 6), np.nan)
    seqs = [None] * npil
    maps = np.zeros((npil, 4), int)
    p = 0
    while p < npil and k < len(I):
        t = I[k]
        if t == -1:
            k += 1
            p += 1
            continue
        if t != 2:
            raise ValueError("unexpected pillar record type %d at pillar %d" % (t, p))
        lines[p] = F[k + 1:k + 7]
        ns = int(I[k + 7])
        k += 8
        if ns > 1:
            maps[p] = I[k:k + 4]
            k += 4
        z = F[k:k + ns * nlev].astype(float).reshape(ns, nlev)
        z[z > FLT_MAX_CUT] = np.nan
        seqs[p] = z
        k += ns * nlev
        p += 1

    # ---- Eclipse COORD / ZCORN (depth positive, metres, Eclipse j = ny-1-r)
    coord = np.full((ny + 1, nx + 1, 6), np.nan)
    for pid in range(npil):
        pr, pi = divmod(pid, nx + 1)
        L = lines[pid]
        if np.isnan(L[0]):
            continue
        coord[ny - pr, pi] = [L[0] + x0, L[1] + y0, -L[2], L[3] + x0, L[4] + y0, -L[5]]
    # fill empty pillars (no geometry) with a vertical line through the nearest neighbour so the
    # array is well-formed; cells touching them are undefined anyway (NaN zcorn)
    if np.isnan(coord).any():
        mask = np.isnan(coord[:, :, 0])
        filled = coord.copy()
        jj, ii = np.nonzero(~mask)
        for j_, i_ in zip(*np.nonzero(mask)):
            n_ = np.argmin((jj - j_) ** 2 + (ii - i_) ** 2)
            filled[j_, i_] = coord[jj[n_], ii[n_]]
        coord = filled

    zc = np.full((2 * nz, 2 * ny, 2 * nx), np.nan)
    # eclipse corner (di,dj) of cell (i,j)  <->  AREAL corner q of GSG column (i, r=ny-1-j)
    q_of = {(0, 0): 0, (1, 0): 3, (0, 1): 1, (1, 1): 2}
    lev = np.arange(nz)
    for r in range(ny):
        j = ny - 1 - r
        for i in range(nx):
            c = r * nx + i
            for (di, dj), qq in q_of.items():
                pid = colp[c, qq]
                z = seqs[pid]
                if z is None:
                    continue
                s = maps[pid][3 - qq] if z.shape[0] > 1 else 0
                depth = -z[s]
                zc[2 * lev, 2 * j + dj, 2 * i + di] = depth[:-1]
                zc[2 * lev + 1, 2 * j + dj, 2 * i + di] = depth[1:]

    props = {}
    if read_props:
        for off in idx.get("PROP", []):
            # body: int32 idx, string name, 'c   ', int32 isfloat, int32 0, int32 count, data
            nm_len = struct.unpack_from("<i", d, off + 4)[0]
            nm = d[off + 8:off + 8 + nm_len].decode("latin1").strip()
            b = off + 8 + nm_len + 4
            isfloat, _, cnt = struct.unpack_from("<iii", d, b)
            if cnt != nx * ny * nz:
                continue
            arr = np.frombuffer(d[b + 12:b + 12 + 4 * cnt], "<f4" if isfloat else "<i4").astype(float)
            arr = arr.reshape(nz, ny, nx)[:, ::-1, :].ravel()  # to Eclipse natural order
            if nm not in props:
                props[nm] = arr
    return dict(name=name, nx=nx, ny=ny, nz=nz, coord=coord.ravel(), zcorn=zc.ravel(), units=units,
                origin=(x0, y0), props=props, n_split_pillars=int(sum(1 for z in seqs if z is not None and z.shape[0] > 1)))

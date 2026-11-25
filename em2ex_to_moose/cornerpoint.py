"""
cornerpoint.py - turn an Eclipse corner-point grid (COORD/ZCORN/ACTNUM) into a conforming-where-
possible HEX8 finite-element mesh, and write it as ExodusII for MOOSE.

Conventions of the produced mesh
  * x, y  : map coordinates (optionally shifted by a constant origin, see Mesh.origin)
  * z     : ELEVATION in metres (z = -depth), so gravity in MOOSE is '0 0 -9.81'
  * nodes are shared between cells that touch the same pillar at the same depth (within tol);
    across faults with throw the two sides do NOT share nodes, i.e. faults are hydraulically
    sealing in a continuous-Galerkin code such as MOOSE PorousFlow (see README).
"""
from __future__ import annotations

import numpy as np

# local corner order (di, dj, dk) for an Exodus HEX8, dk=1 is the deeper (= lower) face
_LOCAL = np.array([(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
                   (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)])
_FLIP = np.array([0, 3, 2, 1, 4, 7, 6, 5])

# Exodus HEX8 side -> logical face, for un-flipped and flipped node ordering
SIDE_OF_FACE = {
    False: {"J-": 1, "I+": 2, "J+": 3, "I-": 4, "K+": 5, "K-": 6},
    True: {"I-": 1, "J+": 2, "I+": 3, "J-": 4, "K+": 5, "K-": 6},
}
# corners (local index, un-flipped order) lying on each logical face
FACE_CORNERS = {"I-": [0, 3, 4, 7], "I+": [1, 2, 5, 6], "J-": [0, 1, 4, 5], "J+": [2, 3, 6, 7],
                "K-": [4, 5, 6, 7], "K+": [0, 1, 2, 3]}
OPPOSITE = {"I-": "I+", "I+": "I-", "J-": "J+", "J+": "J-", "K-": "K+", "K+": "K-"}
STEP = {"I-": (-1, 0, 0), "I+": (1, 0, 0), "J-": (0, -1, 0), "J+": (0, 1, 0),
        "K-": (0, 0, -1), "K+": (0, 0, 1)}


def _hex_dN(xi, eta, zeta):
    """Derivatives of trilinear shape functions at (xi,eta,zeta) in [-1,1]^3 -> (8,3)."""
    s = np.array([(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
                  (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)], float)
    d = np.empty((8, 3))
    d[:, 0] = s[:, 0] * (1 + s[:, 1] * eta) * (1 + s[:, 2] * zeta) / 8
    d[:, 1] = s[:, 1] * (1 + s[:, 0] * xi) * (1 + s[:, 2] * zeta) / 8
    d[:, 2] = s[:, 2] * (1 + s[:, 0] * xi) * (1 + s[:, 1] * eta) / 8
    return d


def hex_jacobians(P, points):
    """P: (n,8,3) node coords. Returns (n, len(points)) determinant of Jacobian."""
    out = np.empty((P.shape[0], len(points)))
    for q, (a, b, c) in enumerate(points):
        dN = _hex_dN(a, b, c)
        J = np.einsum("nki,kj->nij", P, dN)  # (n,3,3): dx_i/dxi_j
        out[:, q] = np.linalg.det(J)
    return out


_G = 1 / np.sqrt(3)
GAUSS8 = [(a, b, c) for a in (-_G, _G) for b in (-_G, _G) for c in (-_G, _G)]


class CornerPointGrid:
    def __init__(self, nx, ny, nz, coord, zcorn, actnum=None, length_scale=1.0):
        self.nx, self.ny, self.nz = int(nx), int(ny), int(nz)
        self.n = self.nx * self.ny * self.nz
        self.coord = np.asarray(coord, float).reshape(self.ny + 1, self.nx + 1, 6) * length_scale
        self.zcorn = np.asarray(zcorn, float).reshape(2 * self.nz, 2 * self.ny, 2 * self.nx) * length_scale
        a = np.ones(self.n, int) if actnum is None else np.asarray(actnum, int)
        self.actnum = a.reshape(self.nz, self.ny, self.nx)

    def gidx(self, i, j, k):
        return (k * self.ny + j) * self.nx + i

    def cell_corners(self, cells):
        """cells: (m,3) int array of (i,j,k), 0-based. Returns (m,8,3) xyz (depth positive) and (m,8) pillar id."""
        i, j, k = cells[:, 0], cells[:, 1], cells[:, 2]
        P = np.empty((len(cells), 8, 3))
        pil = np.empty((len(cells), 8), np.int64)
        for c, (di, dj, dk) in enumerate(_LOCAL):
            z = self.zcorn[2 * k + dk, 2 * j + dj, 2 * i + di]
            pi_, pj = i + di, j + dj
            line = self.coord[pj, pi_]
            xt, yt, zt, xb, yb, zb = line.T
            dz = zb - zt
            with np.errstate(invalid="ignore", divide="ignore"):
                t = np.where(np.abs(dz) > 1e-12, (z - zt) / np.where(dz == 0, 1, dz), 0.0)
            P[:, c, 0] = xt + (xb - xt) * t
            P[:, c, 1] = yt + (yb - yt) * t
            P[:, c, 2] = z
            pil[:, c] = pj * (self.nx + 1) + pi_
        return P, pil


class Mesh:
    """Result of CornerPointGrid -> HEX8 conversion."""

    def __init__(self):
        self.nodes = None          # (N,3)
        self.conn = None           # (E,8) 0-based node ids, Exodus order
        self.cells = None          # (E,3) i,j,k 0-based of each element
        self.gidx = None           # (E,) global (natural) cell index
        self.flipped = None        # (E,) bool
        self.sidesets = {}         # name -> (elem_ids (0-based), side_ids (1-based))
        self.blocks = None         # (E,) int block id
        self.block_names = {}      # id -> name
        self.origin = np.zeros(3)
        self.stats = {}
        self.elem_vars = {}        # name -> (E,) float
        self.node_vars = {}        # name -> (N,) float
        self.info = []


def snap_fault_nodes(nodes, conn, pil, snap=1.0, frac=0.3, region=None, log=print):
    """Merge mesh nodes that sit on the same pillar at almost the same depth on the two sides of a fault.

    A fault with a throw much smaller than the cells it cuts (a few tens of cm on metre-thick layers) makes
    the faces on its two sides slightly offset, so they share no nodes and MOOSE would see no flow across
    them, while Eclipse connects such cells. Nodes of one pillar are merged when their elevations differ by
    less than min(snap, frac x thickness of the thinnest cell using them); the top and bottom of one cell
    are never merged. Faults with larger throw stay disconnected.
    nodes (N,3), conn (E,8), pil (E,8) pillar id of each element corner; returns (nodes, conn, n_merged)."""
    E = len(conn)
    if snap <= 0 or E == 0:
        return nodes, conn, 0
    # thickness of the element along each of its 4 pillars (corners c and c^4 lie on the same pillar)
    z = nodes[conn][:, :, 2]
    pair = np.array([4, 5, 6, 7, 0, 1, 2, 3])
    h = np.abs(z - z[:, pair])
    tol_inst = np.minimum(snap, frac * h)
    N = len(nodes)
    tol = np.full(N, np.inf)
    np.minimum.at(tol, conn.ravel(), tol_inst.ravel())
    npil = -np.ones(N, np.int64)
    npil[conn.ravel()] = pil.ravel()
    if region is not None:
        nreg = np.zeros(N, np.int64)
        nreg[conn.ravel()] = np.repeat(np.asarray(region, np.int64), 8)
    else:
        nreg = np.zeros(N, np.int64)
    used = npil >= 0

    def cluster(block):
        ids = np.nonzero(used & ~block[np.maximum(npil, 0)])[0]
        o = np.lexsort((nodes[ids, 2], nreg[ids], npil[ids]))
        ids = ids[o]
        lab = np.arange(N)
        if len(ids) > 1:
            same = (npil[ids[1:]] == npil[ids[:-1]]) & (nreg[ids[1:]] == nreg[ids[:-1]])
            close = np.diff(nodes[ids, 2]) < np.minimum(tol[ids[1:]], tol[ids[:-1]])
            join = same & close
            start = np.r_[True, ~join]
            cid = np.cumsum(start) - 1
            first = ids[start]
            lab[ids] = first[cid]
        return lab

    block = np.zeros(int(npil.max()) + 1, bool)
    for _ in range(5):
        lab = cluster(block)
        c2 = lab[conn]
        bad = (c2 == c2[:, pair]).any(1)       # an element whose top and bottom corner would merge
        if not bad.any():
            break
        bp = pil[bad][c2[bad] == c2[bad][:, pair]]
        block[np.unique(bp)] = True
    lab = cluster(block)
    uniq, inv = np.unique(lab, return_inverse=True)
    newn = np.zeros((len(uniq), 3))
    cnt = np.zeros(len(uniq))
    np.add.at(newn, inv, nodes)
    np.add.at(cnt, inv, 1)
    keepn = np.isin(np.arange(N), np.unique(conn))
    newn /= np.maximum(cnt, 1)[:, None]
    conn2 = inv[conn]
    # drop nodes no longer referenced
    usedn = np.unique(conn2)
    remap = -np.ones(len(newn), np.int64)
    remap[usedn] = np.arange(len(usedn))
    n_merged = int(len(np.unique(conn)) - len(usedn))
    log("  fault snapping: %d nodes merged across small fault throws (tolerance min(%g m, %g x cell thickness))"
        % (n_merged, snap, frac))
    return newn[usedn], remap[conn2], n_merged


def build_mesh(grid: CornerPointGrid, keep=None, merge_tol=1e-3, min_rel_jac=1e-6,
               shift_origin=True, split_region=None, log=print, fault_snap=0.0, fault_snap_frac=0.3):
    """Convert the grid. keep: optional boolean (n,) mask (in addition to ACTNUM).
    split_region: optional integer (n,) array; cells of different value never share nodes
    (used for zero-transmissibility MULTREGT barriers)."""
    act = grid.actnum.ravel() > 0
    if keep is not None:
        act &= np.asarray(keep, bool)
    gid = np.nonzero(act)[0]
    k, rem = np.divmod(gid, grid.nx * grid.ny)
    j, i = np.divmod(rem, grid.nx)
    cells = np.stack([i, j, k], 1)
    P, pil = grid.cell_corners(cells)
    P[:, :, 2] *= -1.0  # depth -> elevation

    finite = np.isfinite(P).all(axis=(1, 2))
    # orientation: flip lateral order where centre Jacobian is negative
    detc = hex_jacobians(P, [(0.0, 0.0, 0.0)])[:, 0]
    flipped = detc < 0
    order = np.where(flipped[:, None], _FLIP[None, :], np.arange(8)[None, :])
    P = np.take_along_axis(P, order[:, :, None], axis=1)
    pil = np.take_along_axis(pil, order, axis=1)

    # degenerate / inverted cells
    dets = hex_jacobians(P, GAUSS8)
    vol = dets.sum(1)  # Gauss weights are 1 for 2x2x2 rule
    ref = np.median(np.abs(vol[finite])) if finite.any() else 1.0
    good = finite & (dets.min(1) > min_rel_jac * ref / 8.0)
    nbad = int((~good).sum())
    log("  cells: %d in grid, %d active, %d dropped as degenerate/pinched/inverted (%.3f%% of active bulk volume)"
        % (grid.n, len(gid), nbad, 100.0 * np.abs(vol[~good & finite]).sum() / max(np.abs(vol[finite]).sum(), 1e-30)))
    P, pil, cells, gid, flipped, vol = P[good], pil[good], cells[good], gid[good], flipped[good], vol[good]

    # merge nodes: same pillar and same elevation (within merge_tol)
    zq = np.round(P[:, :, 2] / merge_tol).astype(np.int64)
    if split_region is not None:
        reg = np.repeat(np.asarray(split_region, np.int64)[gid], 8)
        key = np.stack([pil.ravel(), zq.ravel(), reg], 1)
    else:
        key = np.stack([pil.ravel(), zq.ravel()], 1)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    conn = inv.reshape(-1, 8)
    nodes = np.zeros((len(uniq), 3))
    cnt = np.zeros(len(uniq))
    np.add.at(nodes, inv, P.reshape(-1, 3))
    np.add.at(cnt, inv, 1)
    nodes /= cnt[:, None]

    # repeated nodes inside an element (partially collapsed) -> drop if any (libMesh needs 8 distinct)
    srt = np.sort(conn, 1)
    distinct = (np.diff(srt, axis=1) != 0).all(1)
    if not distinct.all():
        nd = int((~distinct).sum())
        log("  %d partially-collapsed cells (coincident corners) dropped" % nd)
        conn, cells, gid, flipped, vol = conn[distinct], cells[distinct], gid[distinct], flipped[distinct], vol[distinct]
        used = np.unique(conn)
        remap = -np.ones(len(nodes), np.int64)
        remap[used] = np.arange(len(used))
        nodes, conn = nodes[used], remap[conn]

    n_snap = 0
    if fault_snap > 0:
        pil_d = pil if distinct.all() else pil[distinct]
        reg_d = None if split_region is None else np.asarray(split_region, np.int64)[gid]
        nodes, conn, n_snap = snap_fault_nodes(nodes, conn, pil_d, fault_snap, fault_snap_frac, reg_d, log)

    m = Mesh()
    if shift_origin:
        m.origin = np.array([np.floor(nodes[:, 0].min() / 1000) * 1000,
                             np.floor(nodes[:, 1].min() / 1000) * 1000, 0.0])
        nodes = nodes - m.origin
    m.stats["fault_snapped_nodes"] = n_snap
    m.nodes, m.conn, m.cells, m.gidx, m.flipped = nodes, conn, cells, gid, flipped
    m.stats.update(n_grid=grid.n, n_active=int(act.sum()), n_elem=len(conn), n_nodes=len(nodes),
                   n_flipped=int(flipped.sum()), bulk_volume=float(vol.sum()))

    _sidesets_and_faults(grid, m, log)
    return m


def _sidesets_and_faults(grid, m, log):
    """Logical boundary sidesets + statistics on non-conforming (faulted) internal faces."""
    nx, ny, nz = grid.nx, grid.ny, grid.nz
    lut = -np.ones(grid.n, np.int64)
    lut[m.gidx] = np.arange(len(m.gidx))
    i, j, k = m.cells.T
    # node ids in UN-flipped local order (so FACE_CORNERS applies)
    conn_local = np.where(m.flipped[:, None], m.conn[:, _FLIP], m.conn)
    names = {"I-": "imin", "I+": "imax", "J-": "jmin", "J+": "jmax", "K-": "top", "K+": "bottom"}
    faces = list(STEP)
    # geometric test: a face is on the mesh boundary iff no other element has the same 4 nodes
    keys = np.concatenate([np.sort(conn_local[:, FACE_CORNERS[f]], 1) for f in faces])
    _, inv, cnt = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
    shared = (cnt[inv.ravel()] > 1).reshape(len(faces), -1)
    m.face_external = {}
    fault_e, fault_s = [], []
    nonconf = 0
    internal = 0
    for fi, face in enumerate(faces):
        di, dj, dk = STEP[face]
        ii, jj, kk = i + di, j + dj, k + dk
        inside = (ii >= 0) & (ii < nx) & (jj >= 0) & (jj < ny) & (kk >= 0) & (kk < nz)
        nb = -np.ones(len(i), np.int64)
        nb[inside] = lut[(kk[inside] * ny + jj[inside]) * nx + ii[inside]]
        geo_ext = ~shared[fi]
        ext = geo_ext & (nb < 0)          # true outer boundary (or edge of an inactive region)
        flt = geo_ext & (nb >= 0)         # logical neighbour exists but mesh is not connected: fault
        m.face_external[face] = geo_ext
        for sel, dst in ((ext, None), (flt, "fault")):
            e = np.nonzero(sel)[0]
            side = np.where(m.flipped[e], SIDE_OF_FACE[True][face], SIDE_OF_FACE[False][face])
            if dst is None:
                m.sidesets[names[face]] = (e, side)
            else:
                fault_e.append(e)
                fault_s.append(side)
        if face in ("I+", "J+", "K+"):
            internal += int((nb >= 0).sum())
            nonconf += int(flt.sum())
    m.sidesets["fault_faces"] = (np.concatenate(fault_e), np.concatenate(fault_s))
    m.stats.update(internal_faces=internal, nonconforming_faces=nonconf)
    log("  %d logical internal faces, %d (%.2f%%) non-conforming (fault throw) -> no-flow in MOOSE"
        % (internal, nonconf, 100.0 * nonconf / max(internal, 1)))


def add_face_sideset(grid, m, name, cell_mask, face):
    """Sideset on logical face `face` of kept cells where cell_mask (global, natural order) is true.
    Only faces that are external to the mesh are kept (MOOSE integrated BCs need a boundary)."""
    sel = np.asarray(cell_mask, bool)[m.gidx] & m.face_external[face]
    e = np.nonzero(sel)[0]
    side = np.where(m.flipped[e], SIDE_OF_FACE[True][face], SIDE_OF_FACE[False][face])
    m.sidesets[name] = (e, side)
    return len(e)


def cell_centres(m):
    return m.nodes[m.conn].mean(1)


def element_to_node(m, values, weights=None):
    """Average element values onto nodes (weighted)."""
    w = np.ones(len(m.conn)) if weights is None else weights
    num = np.zeros(len(m.nodes))
    den = np.zeros(len(m.nodes))
    for c in range(8):
        np.add.at(num, m.conn[:, c], values * w)
        np.add.at(den, m.conn[:, c], w)
    return num / np.maximum(den, 1e-300)


# ----------------------------------------------------------------------------------------------
# ExodusII writer (netCDF, 64-bit offset classic format, as written by SEACAS "large model")
# ----------------------------------------------------------------------------------------------
def write_exodus(path, m: Mesh, title="petrel2moose"):
    import netCDF4

    E, N = m.conn.shape[0], m.nodes.shape[0]
    blocks = np.zeros(E, int) if m.blocks is None else np.asarray(m.blocks, int)
    bids = sorted(set(blocks.tolist()))
    # Exodus wants elements numbered contiguously per block
    perm = np.argsort(blocks, kind="stable")
    inv = np.empty(E, np.int64)
    inv[perm] = np.arange(E)
    conn = m.conn[perm]
    bl = blocks[perm]

    ds = netCDF4.Dataset(path, "w", format="NETCDF3_64BIT_OFFSET")
    ds.api_version = np.float32(8.03)
    ds.version = np.float32(8.03)
    ds.floating_point_word_size = np.int32(8)
    ds.file_size = np.int32(1)
    ds.maximum_name_length = np.int32(32)
    ds.int64_status = np.int32(0)
    ds.title = title

    ds.createDimension("len_string", 33)
    ds.createDimension("len_line", 81)
    ds.createDimension("four", 4)
    ds.createDimension("len_name", 33)
    ds.createDimension("time_step", None)
    ds.createDimension("num_dim", 3)
    ds.createDimension("num_nodes", N)
    ds.createDimension("num_elem", E)
    ds.createDimension("num_el_blk", len(bids))

    def names_var(vname, dim, names):
        v = ds.createVariable(vname, "S1", (dim, "len_name"))
        arr = np.zeros((len(names), 33), "S1")
        for r, s in enumerate(names):
            b = s.encode()[:32]
            arr[r, :len(b)] = np.frombuffer(b, dtype="S1")
        v[:] = arr

    if m.info:
        ds.createDimension("num_info", len(m.info))
        v = ds.createVariable("info_records", "S1", ("num_info", "len_line"))
        arr = np.zeros((len(m.info), 81), "S1")
        for r, s in enumerate(m.info):
            b = s.encode()[:80]
            arr[r, :len(b)] = np.frombuffer(b, dtype="S1")
        v[:] = arr

    t = ds.createVariable("time_whole", "f8", ("time_step",))
    t[0] = 0.0

    for c, nm in enumerate("xyz"):
        v = ds.createVariable("coord" + nm, "f8", ("num_nodes",))
        v[:] = m.nodes[:, c]
    names_var("coor_names", "num_dim", ["x", "y", "z"])

    st = ds.createVariable("eb_status", "i4", ("num_el_blk",))
    st[:] = 1
    pr = ds.createVariable("eb_prop1", "i4", ("num_el_blk",))
    pr.setncattr("name", "ID")
    pr[:] = bids
    names_var("eb_names", "num_el_blk", [m.block_names.get(b, "block_%d" % b) for b in bids])
    for n, b in enumerate(bids, 1):
        sel = bl == b
        ds.createDimension("num_el_in_blk%d" % n, int(sel.sum()))
        ds.createDimension("num_nod_per_el%d" % n, 8)
        v = ds.createVariable("connect%d" % n, "i4", ("num_el_in_blk%d" % n, "num_nod_per_el%d" % n))
        v.elem_type = "HEX8"
        v[:] = conn[sel] + 1

    ss = [(nm, e, s) for nm, (e, s) in m.sidesets.items() if len(e)]
    if ss:
        ds.createDimension("num_side_sets", len(ss))
        v = ds.createVariable("ss_status", "i4", ("num_side_sets",))
        v[:] = 1
        v = ds.createVariable("ss_prop1", "i4", ("num_side_sets",))
        v.setncattr("name", "ID")
        v[:] = np.arange(1, len(ss) + 1)
        names_var("ss_names", "num_side_sets", [s[0] for s in ss])
        for n, (nm, e, s) in enumerate(ss, 1):
            ds.createDimension("num_side_ss%d" % n, len(e))
            v = ds.createVariable("elem_ss%d" % n, "i4", ("num_side_ss%d" % n,))
            v[:] = inv[e] + 1
            v = ds.createVariable("side_ss%d" % n, "i4", ("num_side_ss%d" % n,))
            v[:] = s

    # (no elem_num_map: libMesh would use it as element ids and they must be contiguous;
    #  the Eclipse cell of each element is in the cell_i/cell_j/cell_k element variables)

    if m.elem_vars:
        evn = list(m.elem_vars)
        ds.createDimension("num_elem_var", len(evn))
        names_var("name_elem_var", "num_elem_var", evn)
        tab = ds.createVariable("elem_var_tab", "i4", ("num_el_blk", "num_elem_var"))
        tab[:] = 1
        for vi, nm in enumerate(evn, 1):
            vals = np.asarray(m.elem_vars[nm], float)[perm]
            for n, b in enumerate(bids, 1):
                v = ds.createVariable("vals_elem_var%deb%d" % (vi, n), "f8", ("time_step", "num_el_in_blk%d" % n))
                v[0, :] = vals[bl == b]
    if m.node_vars:
        nvn = list(m.node_vars)
        ds.createDimension("num_nod_var", len(nvn))
        names_var("name_nod_var", "num_nod_var", nvn)
        for vi, nm in enumerate(nvn, 1):
            v = ds.createVariable("vals_nod_var%d" % vi, "f8", ("time_step", "num_nodes"))
            v[0, :] = m.node_vars[nm]
    ds.close()
    return perm

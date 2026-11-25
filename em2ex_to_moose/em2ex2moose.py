#!/usr/bin/env python3
"""em2ex2moose.py - build a MOOSE PorousFlow input file for an Eclipse model on a mesh made by em2ex.

The grid (nodes, elements) and the cell properties (porosity, NTG, permeability, SATNUM) come from
the em2ex Exodus file. Everything else comes from the Eclipse case, using the same engine as
petrel2moose.py:
  * gas and water PVT (PVDG, PVTW, DENSITY), relative permeability (one Corey fit per SATNUM table),
  * initial state: saturations from the restart file (.X0000 / .UNRST), pressure hydrostatic from EQUIL,
  * Fetkovich / Carter-Tracy aquifers (AQUANCON faces become sidesets),
  * well connections (COMPDAT) and the rate history (WCONHIST / WCONPROD / ...).

It writes:
  <out>.e              the em2ex mesh (same nodes and elements) plus initial-state fields and sidesets
  <out>.i              the MOOSE input file
  <out>_gas_pvt.csv    gas PVT table read by the input file
  <out>_wells.csv      well connections
  <out>_report.md      every step and approximation
  <out>_meta.json      origin shift, block names, sideset sizes

Example (Croft; em2ex run first with the same --coarsen-z):
  python em2ex.py -f --convert-to-m --coarsen-z 5 --coarsen-z-pinch -o croft_19.e croft.grdecl
  python em2ex2moose.py --mesh croft_19.e --case CROFT_files/HM_MODEL.DATA --coarsen-z 5 -o croft_em2ex

  then:  mpiexec -n 2 porous_flow-opt -i croft_em2ex.i

The em2ex mesh must be written in metres (em2ex --convert-to-m for FIELD grids). Depth may be
positive (em2ex default) or negative (em2ex --flip). em2ex's --translate, --mapaxes and --refine-xy
change the geometry and are not supported here.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import netCDF4  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

import cornerpoint as cp  # noqa: E402
import petrel2moose as p2m  # noqa: E402

RB_TO_M3 = 0.158987294928


# ----------------------------------------------------------------------------------------------
# reading the em2ex mesh
# ----------------------------------------------------------------------------------------------
def read_em2ex(path):
    """Nodes (N,3), connectivity (E,8, 0-based) and element variables of an em2ex Exodus file."""
    d = netCDF4.Dataset(path)
    d.set_auto_mask(False)
    X = np.c_[d["coordx"][:], d["coordy"][:], d["coordz"][:]].astype(float)
    nblk = d.dimensions["num_el_blk"].size
    conn = [d["connect%d" % (b + 1)][:].astype(np.int64) - 1 for b in range(nblk)]
    C = np.vstack(conn)
    names = []
    if "name_elem_var" in d.variables:
        names = ["".join(ch.decode() for ch in r if ch) for r in d["name_elem_var"][:]]
    ev = {}
    for k, nm in enumerate(names):
        ev[nm.lower()] = np.concatenate([d["vals_elem_var%deb%d" % (k + 1, b + 1)][0] for b in range(nblk)])
    d.close()
    return X, C, ev


def parse_groups(spec, nz):
    """--coarsen-z value (N or N,N,...) -> list of layer counts that adds up to nz (as em2ex does)."""
    if spec in (None, "", "0", "1"):
        return None
    sizes = [int(s) for s in str(spec).split(",") if s.strip()]
    if len(sizes) == 1:
        n = sizes[0]
        groups = [n] * (nz // n) + ([nz % n] if nz % n else [])
    else:
        groups = sizes
    if sum(groups) != nz:
        raise SystemExit("--coarsen-z %s adds up to %d layers, but the grid has %d" % (spec, sum(groups), nz))
    return groups


# ----------------------------------------------------------------------------------------------
# vertical merging with arbitrary group sizes (same rules as petrel2moose.upscale_layers), used for
# the restart saturations/pressure, SWL and region arrays not carried by the em2ex mesh, and to get
# the merged grid that the em2ex elements are matched against.
# ----------------------------------------------------------------------------------------------
def upscale_groups(md, groups, log=print):
    g = md.grid
    nx, ny, nz = g.nx, g.ny, g.nz
    starts = np.concatenate(([0], np.cumsum(groups)[:-1])).astype(int)
    grp = [(int(s), int(s + n - 1)) for s, n in zip(starts, groups)]
    nzc = len(grp)
    zc = g.zcorn
    zcc = np.empty((2 * nzc, 2 * ny, 2 * nx))
    for kc, (k0, k1) in enumerate(grp):      # top of first layer, bottom of last (as em2ex)
        zcc[2 * kc] = zc[2 * k0]
        zcc[2 * kc + 1] = zc[2 * k1 + 1]
    act = g.actnum.reshape(nz, ny, nx) > 0
    P_ = md.props
    gid = np.nonzero(act.ravel())[0]
    k, rem = np.divmod(gid, nx * ny)
    j, i = np.divmod(rem, nx)
    Pc, _ = g.cell_corners(np.stack([i, j, k], 1))
    vol = np.zeros(g.n)
    vol[gid] = np.abs(p2m._hex_volume(Pc))
    porv = np.nan_to_num(P_["PORV"]) if "PORV" in P_ else \
        np.nan_to_num(P_["PORO"]) * np.nan_to_num(P_.get("NTG", np.ones(g.n)), nan=1.0) * vol
    porv = np.where(act.ravel(), porv, 0.0)

    def R(a):
        return np.asarray(a, float).reshape(nz, ny, nx)

    def agg(fn):
        return np.stack([fn(slice(k0, k1 + 1)) for k0, k1 in grp]).ravel()

    A, PV = act, R(porv)
    pvsum = agg(lambda s_: PV[s_].sum(0))
    newp = {"PORV": pvsum}
    for nm in ("SWL", "SWATINIT", "PCW", "SWLPC", "KRW", "KRG", "IKRW", "IKRG"):
        if nm in P_:
            x = R(np.nan_to_num(P_[nm]))
            newp[nm] = agg(lambda s_: (x[s_] * PV[s_]).sum(0)) / np.maximum(pvsum, 1e-30)
    for nm in ("SATNUM", "EQLNUM", "FIPNUM", "PVTNUM", "MULTNUM", "FLUXNUM", "OPERNUM", "IMBNUM"):
        if nm in P_:
            x = R(np.nan_to_num(P_[nm])).astype(int)
            vals = np.unique(x[A])

            def maj(s_, x=x, vals=vals):
                best = np.zeros(x.shape[1:])
                arg = np.full(x.shape[1:], vals[0])
                for v in vals:
                    w_ = (PV[s_] * (x[s_] == v) * A[s_]).sum(0)
                    arg = np.where(w_ > best, v, arg)
                    best = np.maximum(best, w_)
                return arg
            newp[nm] = agg(maj).astype(float)
    actc = agg(lambda s_: A[s_].any(0)).astype(int)
    rst = {}
    for nm, val in (md.restart or {}).items():
        x = R(np.nan_to_num(val))
        rst[nm] = agg(lambda s_: (x[s_] * PV[s_]).sum(0)) / np.maximum(pvsum, 1e-30)
    out = p2m.Model()
    out.units, out.deck, out.mapaxes = md.units, md.deck, md.mapaxes
    out.notes = list(md.notes)
    out.grid = p2m.CornerPointGrid(nx, ny, nzc, g.coord.ravel(), zcc.ravel(), actc)
    out.props, out.restart = newp, rst
    kmap = np.zeros(nz + 1, int)
    for kc, (k0, k1) in enumerate(grp):
        kmap[k0 + 1:k1 + 2] = kc + 1
    out.kmap = kmap
    out.fine_porv_total = float(porv.sum())
    log("  layers merged as in em2ex: %d -> %d (group sizes %s)" % (nz, nzc, " ".join(map(str, groups))))
    return out


# ----------------------------------------------------------------------------------------------
# em2ex elements -> cells of the (merged) Eclipse grid, in petrel2moose's mesh layout
# ----------------------------------------------------------------------------------------------
def mesh_from_em2ex(path, grid, shift_origin=True, log=print, fault_snap=0.0, fault_snap_frac=0.3):
    X, C, ev = read_em2ex(path)
    E = len(C)
    # depth sign: em2ex writes z = depth (positive) unless --flip was used
    zsign = 1.0 if np.median(X[:, 2]) > 0 else -1.0
    Xd = X.copy()
    Xd[:, 2] *= zsign                                  # now depth positive, as in the Eclipse grid

    act = grid.actnum.ravel() > 0
    gid = np.nonzero(act)[0]
    k, rem = np.divmod(gid, grid.nx * grid.ny)
    j, i = np.divmod(rem, grid.nx)
    cells = np.stack([i, j, k], 1)
    P, pil = grid.cell_corners(cells)                  # (n,8,3) depth positive, _LOCAL order
    cen_g = P.mean(1)
    cen_e = Xd[C].mean(1)
    span = np.ptp(cen_g, axis=0).max()
    d, idx = cKDTree(cen_g).query(cen_e)
    bad = d > 0.5                                      # metres
    if bad.any() or len(np.unique(idx)) != E:
        ratio = np.ptp(cen_e, axis=0).max() / max(span, 1e-30)
        hint = ""
        if abs(ratio - 3.2808) < 0.1:
            hint = " The em2ex mesh looks like it is in feet: rerun em2ex with --convert-to-m."
        raise SystemExit(
            "The em2ex mesh does not line up with the Eclipse grid (%d of %d elements more than 0.5 m from any "
            "cell centre, worst %.1f m).%s Check that --coarsen-z here is the same as in em2ex, and that em2ex "
            "was run without --translate / --mapaxes / --refine-xy / --extract-*." % (int(bad.sum()), E, d.max(), hint))
    log("  em2ex mesh: %d elements matched to Eclipse cells (largest centre offset %.2g m)" % (E, d.max()))

    cells_e, gid_e, P_e = cells[idx], gid[idx], P[idx]
    # put each element's nodes in petrel2moose's local corner order: node of the em2ex element
    # nearest to each true corner of the matched cell
    loc = np.empty((E, 8), np.int64)
    for c in range(8):
        dist = np.linalg.norm(Xd[C] - P_e[:, c:c + 1, :], axis=2)   # (E,8)
        loc[:, c] = np.argmin(dist, 1)
    ok = np.array([len(set(r)) == 8 for r in loc])
    if not ok.all():
        raise SystemExit("%d em2ex elements have corners that cannot be matched to their Eclipse cell" % int((~ok).sum()))
    conn = np.take_along_axis(C, loc, 1)
    nodes = Xd.copy()
    nodes[:, 2] *= -1.0                                # elevation (z = -depth), as petrel2moose
    Pl = nodes[conn]
    detc = cp.hex_jacobians(Pl, [(0.0, 0.0, 0.0)])[:, 0]
    flipped = detc < 0
    order = np.where(flipped[:, None], cp._FLIP[None, :], np.arange(8)[None, :])
    conn = np.take_along_axis(conn, order, 1)

    pil_e = np.take_along_axis(pil[idx], order, 1)     # pillar of each corner, in conn order
    n_snap = 0
    if fault_snap > 0:
        nodes, conn, n_snap = cp.snap_fault_nodes(nodes, conn, pil_e, fault_snap, fault_snap_frac, None, log)
    m = cp.Mesh()
    m.stats["fault_snapped_nodes"] = n_snap
    if shift_origin:
        m.origin = np.array([np.floor(nodes[:, 0].min() / 1000) * 1000,
                             np.floor(nodes[:, 1].min() / 1000) * 1000, 0.0])
        nodes = nodes - m.origin
    m.nodes, m.conn, m.cells, m.gidx, m.flipped = nodes, conn, cells_e, gid_e, flipped
    vol = np.abs(p2m._hex_volume(nodes[conn]))
    m.stats.update(n_grid=grid.n, n_active=int(act.sum()), n_elem=E, n_nodes=len(nodes),
                   n_flipped=int(flipped.sum()), bulk_volume=float(vol.sum()))
    cp._sidesets_and_faults(grid, m, log)
    m.em2ex = dict(path=path, ev=ev, zsign=zsign, n_missing=int(act.sum()) - E)
    return m


def em2ex_props(m, n, ev):
    """Global (natural order) property arrays taken from the em2ex element variables."""
    out = {}
    lookup = {"PORO": "poro", "NTG": "ntg", "PERMX": "permx", "PERMY": "permy", "PERMZ": "permz",
              "SATNUM": "satnum", "FIPNUM": "fipnum", "PVTNUM": "pvtnum", "EQLNUM": "eqlnum"}
    for key, nm in lookup.items():
        if nm in ev:
            a = np.zeros(n)
            a[m.gidx] = ev[nm]
            out[key] = a
    return out


# ----------------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="MOOSE PorousFlow input for an Eclipse case on a mesh made by em2ex",
        epilog="Other petrel2moose options can be added, e.g. --restart-pressure, --mobile-connate-water, "
               "--temperature, --no-shift.")
    ap.add_argument("--mesh", required=True, help="Exodus file written by em2ex (in metres)")
    ap.add_argument("--case", required=True, help="Eclipse .DATA deck (EGRID, INIT and restart next to it)")
    ap.add_argument("--coarsen-z", default=None,
                    help="the same --coarsen-z value that was given to em2ex (N or N,N,...); omit if none")
    ap.add_argument("-o", "--output", default="model", help="output prefix (default: model)")
    ap.add_argument("--restart-step", type=int, default=0,
                    help="restart step for the initial saturations (default 0); -1 = no restart, use EQUIL")
    ap.add_argument("--average-rates", type=float, default=30.0,
                    help="average well rates over periods of this many days (default 30; 0 = exact schedule)")
    ap.add_argument("--capillary", choices=("none", "deck"), default="none",
                    help="none: Pc = 0 with immobile connate water (default, the validated setup); deck: capillary "
                         "pressure from the deck's Pc tables, scaled cell by cell like Eclipse (needed for fine grids)")
    ap.add_argument("--fv", action="store_true",
                    help="cell-centred finite-volume input (one saturation per cell, like Eclipse) instead of the "
                         "default nodal finite-element one; uses the FVPorousFlowAdvectiveFluxTPFA add-on kernel")
    ap.add_argument("--older-moose", action="store_true",
                    help="leave out time_integration_scheme (MOOSE builds before late 2026 reject it and "
                         "integrate the aquifer influx with the trapezoidal rule instead)")
    ap.add_argument("--no-checkpoint", action="store_true",
                    help="do not write MOOSE checkpoints (by default they are written, so a stopped run can "
                         "be continued with --recover)")
    args, rest = ap.parse_known_args(argv)
    state = {}

    # How petrel2moose.convert() runs, and where em2ex comes in:
    #   md = load_from_egrid(...)            -> load the Eclipse case, merge layers as em2ex did,
    #                                           match the em2ex elements to the merged cells and
    #                                           replace PORO/NTG/PERM*/SATNUM by em2ex's values
    #   md = upscale_layers(md, ...)         -> hand back that merged model (keeps the layer map
    #                                           used for wells and aquifer connections)
    #   m = build_mesh(grid, ...)            -> the em2ex mesh instead of petrel2moose's own
    #   ... PVT, rel-perm, initial state, aquifers, wells, mesh and .i writing: unchanged
    orig_load = p2m.load_from_egrid

    def load(*a, **kw):
        md = orig_load(*a, **kw)
        nz_fine = md.grid.nz
        fine_act = md.grid.actnum.ravel() > 0
        fine_pv = float(np.nan_to_num(md.props.get("PORV", np.zeros(md.grid.n)))[fine_act].sum())
        state["fine_pv_m3"] = fine_pv * (RB_TO_M3 if md.units == "FIELD" else 1.0)
        groups = parse_groups(args.coarsen_z, nz_fine)
        state.update(fine_n_layers=nz_fine, groups=groups)
        if groups:
            md = upscale_groups(md, groups, print)
        m = mesh_from_em2ex(args.mesh, md.grid, log=print)
        props = em2ex_props(m, md.grid.n, m.em2ex["ev"])
        if "PORO" not in props or "PERMX" not in props:
            raise SystemExit("The em2ex mesh has no poro / permx element variables")
        inmesh = np.zeros(md.grid.n, bool)
        inmesh[m.gidx] = True
        md.props.update(props)
        md.props.pop("PORV", None)             # pore volume = em2ex PORO x NTG x cell volume
        md.props["PORO"] = np.where(inmesh, md.props["PORO"], 0.0)
        state.update(md=md, kmap=getattr(md, "kmap", None), em2ex_missing=m.em2ex["n_missing"])
        return md

    def upscale_layers(md, nmerge, log=print):
        md = state["md"]
        md.kmap = state["kmap"]
        return md

    def build_mesh(grid, keep=None, merge_tol=1e-3, shift_origin=True, split_region=None, log=print,
                   fault_snap=0.0, fault_snap_frac=0.3):
        if keep is not None:
            n_drop = int(np.sum(~np.asarray(keep, bool)[grid.actnum.ravel() > 0]))
            state["n_keep_drop"] = n_drop
        return mesh_from_em2ex(args.mesh, grid, shift_origin=shift_origin, log=lambda *x: None,
                               fault_snap=fault_snap, fault_snap_frac=fault_snap_frac)

    orig_convert = p2m.convert

    def convert(a, log=print):
        res = orig_convert(a, log)
        state["params"] = res[1]
        return res

    p2m.load_from_egrid = load
    p2m.upscale_layers = upscale_layers
    p2m.build_mesh = build_mesh
    p2m.convert = convert

    p2m_argv = ["--case", args.case, "-o", args.output, "--average-rates", str(args.average_rates),
                "--capillary", args.capillary] + (["--fv"] if args.fv else [])
    if args.restart_step is not None and args.restart_step >= 0:
        p2m_argv += ["--restart-step", str(args.restart_step)]
    if args.coarsen_z not in (None, "", "0", "1"):
        p2m_argv += ["--merge-layers", "2"]   # only switches on the layer remapping; groups come from --coarsen-z
    p2m.main(p2m_argv + rest)
    _finish_report(args, state)
    _add_meta(args, state)
    if not args.no_checkpoint:
        _add_checkpoint(args.output + ".i")
    if args.older_moose:
        _drop_time_integration_scheme(args.output + ".i")
    _unix_line_endings(args.output)
    print("\nWritten: %s.i  %s.e  %s_gas_pvt.csv  %s_wells.csv  %s_report.md" % ((args.output,) * 5))
    print("Run:     mpiexec -n 2 porous_flow-opt -i %s.i" % os.path.basename(args.output))


def _unix_line_endings(prefix):
    """Write the text outputs with LF line endings, also on Windows: MOOSE on Linux rejects CRLF files."""
    for suffix in (".i", "_gas_pvt.csv", "_wells.csv", "_report.md", "_meta.json"):
        path = prefix + suffix
        if os.path.exists(path):
            with open(path, "rb") as f:
                data = f.read()
            if b"\r\n" in data:
                with open(path, "wb") as f:
                    f.write(data.replace(b"\r\n", b"\n"))


def _add_meta(args, state):
    """Record what the comparison script needs: unit system, surface densities, wells, source files."""
    mp = args.output + "_meta.json"
    if not os.path.exists(mp) or "params" not in state:
        return
    meta = json.load(open(mp))
    prm = state["params"]
    md = state.get("md")
    deck = md.deck if md is not None else None
    rho_o, rho_w_sc, rho_g_sc = deck.densities() if deck else (None, prm.get("rho_w_sc"), None)
    meta.update(units_in=getattr(md, "units", meta.get("units_in")), rho_gas_sc=rho_g_sc, rho_water_sc=rho_w_sc,
                wells=[w["name"] for w in prm.get("wells", [])], em2ex_mesh=os.path.abspath(args.mesh),
                eclipse_case=os.path.abspath(args.case), coarsen_z=args.coarsen_z)
    json.dump(meta, open(mp, "w"), indent=2)


def _add_checkpoint(path):
    """Turn on MOOSE checkpoints in [Outputs] (continue a stopped run with: porous_flow-opt -i X.i --recover)."""
    txt = open(path).read()
    if "checkpoint" in txt or "[Outputs]" not in txt:
        return
    open(path, "w").write(txt.replace("[Outputs]\n", "[Outputs]\n  checkpoint = true\n", 1))


def _drop_time_integration_scheme(path):
    """Remove time_integration_scheme lines, which older MOOSE builds reject as unused parameters."""
    lines = open(path).read().splitlines(True)
    open(path, "w").write("".join(l for l in lines if "time_integration_scheme" not in l))


def _finish_report(args, state):
    """Reword the petrel2moose report lines that describe the mesh and layer merging."""
    rp = args.output + "_report.md"
    if not os.path.exists(rp):
        return
    out = []
    for ln in open(rp).read().splitlines():
        if ln.startswith("# petrel2moose conversion report"):
            ln = "# em2ex2moose conversion report"
        elif ln.startswith("- Vertical upscaling:"):
            g = state.get("groups") or []
            ln = ("- Layers merged exactly as in the em2ex mesh (--coarsen-z %s: %d -> %d layers). Restart "
                  "saturations, SWL and region arrays not in the mesh are pore-volume weighted over the same "
                  "groups; wells and aquifer connections remapped." % (args.coarsen_z, state.get("fine_n_layers", 0), len(g)))
        elif ln.startswith("- Porosity passed to MOOSE is PORO x NTG"):
            ln = "- Porosity passed to MOOSE is em2ex's PORO x NTG (Eclipse pore-volume convention)."
        elif ln.startswith("- Mesh:"):
            out.append("- Mesh geometry from em2ex: %s (same nodes and elements; z converted to elevation, "
                       "x/y shifted by the origin below)." % os.path.basename(args.mesh))
            out.append("- Porosity, NTG, PERMX/Y/Z and SATNUM from the em2ex mesh's element variables.")
        elif ln.startswith("- Pore volume in MOOSE mesh:") and state.get("fine_pv_m3"):
            import re
            mm = re.search(r"([0-9.eE+-]+) m3", ln)
            if mm:
                pv = float(mm.group(1))
                ln = "- Pore volume in MOOSE mesh: %.6g m3 = %.4f%% of Eclipse PORV (%.6g m3)." % (
                    pv, 100 * pv / state["fine_pv_m3"], state["fine_pv_m3"])
        out.append(ln)
    if state.get("em2ex_missing"):
        out.append("- NOTE: %d active cells of the merged Eclipse grid are not in the em2ex mesh and are left out."
                   % state["em2ex_missing"])
    if state.get("n_keep_drop"):
        out.append("- NOTE: %d em2ex elements have PORO x NTG below --min-poro; they stay in the mesh."
                   % state["n_keep_drop"])
    open(rp, "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()

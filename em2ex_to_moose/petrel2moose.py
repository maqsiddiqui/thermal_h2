#!/usr/bin/env python3
"""
petrel2moose.py - convert a Petrel dynamic model (exported/ran as an Eclipse case) into
                  a MOOSE PorousFlow model: ExodusII mesh with cell properties + input file.

A Petrel .pet file is only the project *index*; the grid geometry and property arrays live in
the companion <project>.ptd folder in a proprietary format.  What MOOSE needs is therefore taken
from the Eclipse files Petrel already writes when it runs a simulation case
(<project>.sim/<CASE>/<CASE>.DATA, .EGRID, .INIT, .UNRST ...) or from a GRDECL export.

Typical use
-----------
    # straight from the simulation case Petrel wrote (recommended)
    python petrel2moose.py --case "Penryn_dynamic_model.sim/HM/HM.DATA" -o penryn

    # from a Petrel "Export > Eclipse grid (GRDECL)" with properties in the same/extra files
    python petrel2moose.py --grdecl Penryn.GRDECL --props PORO.GRDECL PERMX.GRDECL -o penryn

Outputs (prefix = -o):
    <prefix>.e          ExodusII mesh: HEX8, z = elevation [m], element vars poro, perm_xx.. [m2],
                        sidesets top/bottom/imin/imax/jmin/jmax (+ one per analytic aquifer)
    <prefix>.i          MOOSE PorousFlow 2-phase (water / methane) input file
    <prefix>_wells.csv  per-connection well table (cell, xyz, rate share)
    <prefix>_report.md  what was converted, what was approximated, and what was dropped

Only numpy + netCDF4 are needed (pip install numpy netCDF4).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from collections import OrderedDict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eclio  # noqa: E402
from cornerpoint import (CornerPointGrid, add_face_sideset, build_mesh, cell_centres,  # noqa: E402
                         element_to_node, write_exodus)

G = 9.81
R_GAS = 8.314462618
M_CH4 = 0.01604246


# ----------------------------------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------------------------------
def _find_case_file(data_path, ext):
    base = os.path.splitext(data_path)[0]
    d = os.path.dirname(os.path.abspath(data_path))
    for cand in (base + ext, base + ext.lower()):
        if os.path.exists(cand):
            return cand
    stem = os.path.basename(base).lower()
    for f in os.listdir(d):
        if f.lower() == stem + ext.lower():
            return os.path.join(d, f)
    return None


class Model:
    """Everything we know about the Eclipse case, in SI units, global (natural) cell order."""

    def __init__(self):
        self.grid = None
        self.props = {}       # name -> (n,) float, global ordering
        self.units = "METRIC"
        self.restart = {}     # PRESSURE [Pa], SGAS, SWAT (global, nan where inactive)
        self.deck = None
        self.notes = []
        self.mapaxes = None


def load_from_egrid(egrid, init=None, unrst=None, restart_step=0, log=print):
    md = Model()
    items = eclio.read_unformatted(egrid, stop_after="ENDGRID")
    gh = eclio.first(items, "GRIDHEAD")
    nx, ny, nz = int(gh[1]), int(gh[2]), int(gh[3])
    unit = eclio.first(items, "GRIDUNIT")
    L = 0.3048 if unit is not None and len(unit) and str(unit[0]).upper().startswith("F") else 1.0
    if unit is not None and len(unit) and str(unit[0]).upper().startswith("C"):
        L = 0.01
    actnum = eclio.first(items, "ACTNUM")
    md.mapaxes = eclio.first(items, "MAPAXES")
    md.grid = CornerPointGrid(nx, ny, nz, eclio.first(items, "COORD"), eclio.first(items, "ZCORN"),
                              actnum, length_scale=L)
    log("  EGRID %s: %d x %d x %d = %d cells, %d active" % (os.path.basename(egrid), nx, ny, nz,
                                                            md.grid.n, int((md.grid.actnum > 0).sum())))
    _attach_init_restart(md, init, unrst, restart_step, log)
    return md


def load_from_gsg(gsg_path, init=None, unrst=None, restart_step=0, log=print):
    """Petrel PETGRID (.GSG) geometry; active cells and properties from the INIT file."""
    import gsg
    G = gsg.read_gsg(gsg_path)
    md = Model()
    nx, ny, nz = G["nx"], G["ny"], G["nz"]
    n = nx * ny * nz
    actnum = None
    if init:
        it = eclio.read_unformatted(init)
        porv = eclio.first(it, "PORV")
        if porv is not None and len(porv) == n:
            actnum = (porv > 0).astype(int)
    if actnum is None:
        a = G["props"].get("ACTNUM_ALL", G["props"].get("ACTNUM"))
        actnum = (np.nan_to_num(a) > 0).astype(int) if a is not None else None
    md.grid = CornerPointGrid(nx, ny, nz, G["coord"], G["zcorn"], actnum, length_scale=1.0)
    # cells without geometry in the GSG cannot be active
    zc = md.grid.zcorn
    undefined = np.isnan(zc).reshape(nz, 2, ny, 2, nx, 2).any(axis=(1, 3, 5)).ravel()
    bad = (md.grid.actnum.ravel() > 0) & undefined
    if bad.any():
        md.notes.append("%d active cells have no geometry in the GSG and were dropped" % int(bad.sum()))
        md.grid.actnum.ravel()[bad] = 0
    md.gsg = G
    log("  GSG %s: grid '%s' %d x %d x %d, %d split (faulted) pillars, %d active cells"
        % (os.path.basename(gsg_path), G["name"], nx, ny, nz, G["n_split_pillars"], int((md.grid.actnum > 0).sum())))
    _attach_init_restart(md, init, unrst, restart_step, log)
    if not init:
        for k, a in G["props"].items():
            md.props.setdefault(k, a)
    # self-check: decoded geometry must reproduce Eclipse's own TOPS / DEPTH / DZ
    if "TOPS" in md.props and "DZ" in md.props:
        L = eclio.UNITS[md.units]["length"]
        act = np.nonzero(md.grid.actnum.ravel() > 0)[0]
        k, rem = np.divmod(act, nx * ny)
        j, i = np.divmod(rem, nx)
        Pc, _ = md.grid.cell_corners(np.stack([i, j, k], 1))
        errs = {}
        for nm, val in (("TOPS", Pc[:, 4:, 2].mean(1)), ("DEPTH", Pc[:, :, 2].mean(1)),
                        ("DZ", (Pc[:, :4, 2] - Pc[:, 4:, 2]).mean(1))):
            if nm in md.props:
                errs[nm] = float(np.nanmax(np.abs(val - md.props[nm][act] * L)))
        md.notes.append("GSG decoding check against Eclipse INIT on %d active cells: max |error| %s"
                        % (len(act), ", ".join("%s %.4f m" % kv for kv in errs.items())))
        log("  GSG self-check vs INIT: " + ", ".join("%s %.4f m" % kv for kv in errs.items()))
        if max(errs.values()) > 0.05:
            log("  WARNING: GSG geometry does not reproduce INIT - this Petrel version may encode GSG differently")
    return md


def _attach_init_restart(md, init, unrst, restart_step, log):
    act = md.grid.actnum.ravel() > 0
    nact = int(act.sum())

    def to_global(a):
        a = np.asarray(a, float)
        if len(a) == md.grid.n:
            return a
        if len(a) == nact:
            g = np.full(md.grid.n, np.nan)
            g[act] = a
            return g
        return None

    if init:
        it = eclio.read_unformatted(init)
        ih = eclio.first(it, "INTEHEAD")
        if ih is not None and len(ih) > 2:
            md.units = {1: "METRIC", 2: "FIELD", 3: "LAB"}.get(int(ih[2]), "METRIC")
        for kw, a in it:
            if kw in ("PORO", "PERMX", "PERMY", "PERMZ", "NTG", "PORV", "SATNUM", "PVTNUM", "FIPNUM",
                      "EQLNUM", "MULTX", "MULTY", "MULTZ", "SWL", "SGCR", "SWCR", "SWU", "SGU", "DEPTH",
                      "DZ", "TOPS", "TRANX", "TRANY", "TRANZ", "ROCKNUM", "MULTPV", "MULTNUM", "FLUXNUM",
                      "OPERNUM", "SWATINIT", "IMBNUM", "PCW", "SWLPC", "KRW", "KRG", "IKRW", "IKRG") and kw not in md.props:
                g = to_global(a)
                if g is not None:
                    md.props[kw] = g * (eclio.UNITS[md.units]["pressure"] if kw == "PCW" else 1.0)
        log("  INIT %s: %s (units %s)" % (os.path.basename(init), " ".join(sorted(md.props)), md.units))
    if unrst:
        md.restart = _read_restart(unrst, restart_step, to_global, md.units, log)


def _read_restart(unrst, step, to_global, units, log):
    items = eclio.read_unformatted(unrst)
    seq_pos = [n for n, (k, _) in enumerate(items) if k == "SEQNUM"]
    if not seq_pos:
        seq_pos = [0]
    if step < 0:
        step = len(seq_pos) + step
    s0 = seq_pos[min(step, len(seq_pos) - 1)]
    s1 = seq_pos[step + 1] if step + 1 < len(seq_pos) else len(items)
    out = {}
    P = eclio.UNITS[units]["pressure"]
    for k, a in items[s0:s1]:
        if k in ("PRESSURE", "SGAS", "SWAT", "PCOW", "PCGW"):
            g = to_global(a)
            if g is not None:
                out[k] = g * (P if k in ("PRESSURE", "PCOW", "PCGW") else 1.0)
    seq = items[s0][1][0] if items[s0][0] == "SEQNUM" else 0
    log("  UNRST report %d (step index %d): %s" % (seq, step, " ".join(out)))
    return out


def load_from_grdecl(grdecl, extra=(), log=print):
    md = Model()
    kws = eclio.read_grdecl(grdecl)
    for e in extra:
        kws.update(eclio.read_grdecl(e))
    if "SPECGRID" in kws:
        nx, ny, nz = [int(v) for v in kws["SPECGRID"][:3]]
    elif "DIMENS" in kws:
        nx, ny, nz = [int(v) for v in kws["DIMENS"][:3]]
    else:
        raise SystemExit("GRDECL has neither SPECGRID nor DIMENS")
    unit = kws.get("GRIDUNIT")
    L = 0.3048 if unit is not None and len(unit) and str(unit[0]).upper().startswith("F") else 1.0
    if L != 1.0:
        md.units = "FIELD"
    md.mapaxes = kws.get("MAPAXES")
    md.grid = CornerPointGrid(nx, ny, nz, kws["COORD"], kws["ZCORN"], kws.get("ACTNUM"), length_scale=L)
    n = md.grid.n
    for k, a in kws.items():
        if k in ("COORD", "ZCORN", "ACTNUM", "SPECGRID", "DIMENS", "MAPAXES", "GRIDUNIT", "MAPUNITS",
                 "COORDSYS"):
            continue
        if a.dtype.kind == "f" and len(a) == n:
            md.props[k] = a
    log("  GRDECL: %d x %d x %d, properties: %s" % (nx, ny, nz, " ".join(md.props)))
    return md


def upscale_layers(md, nmerge, log=print):
    """Merge every `nmerge` consecutive layers into one (vertical upscaling).

    Pore volume (PORV), restart saturations/pressure are pore-volume conserved; horizontal permeability
    is thickness-weighted arithmetic (kh), vertical permeability thickness-weighted harmonic; region
    numbers take the pore-volume majority. Returns the coarse Model; md.kmap maps fine k (1-based) to
    coarse k (1-based)."""
    g = md.grid
    nx, ny, nz = g.nx, g.ny, g.nz
    groups = [(k0, min(k0 + nmerge, nz) - 1) for k0 in range(0, nz, nmerge)]
    nzc = len(groups)
    zc = g.zcorn
    zcc = np.full((2 * nzc, 2 * ny, 2 * nx), np.nan)
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            for kc, (k0, k1) in enumerate(groups):
                zcc[2 * kc] = np.nanmin(zc[2 * k0:2 * k1 + 1:2], axis=0)
                zcc[2 * kc + 1] = np.nanmax(zc[2 * k0 + 1:2 * k1 + 2:2], axis=0)
    act = g.actnum.reshape(nz, ny, nx) > 0
    # fine bulk volumes and thickness for weighting
    gid = np.nonzero(act.ravel())[0]
    k, rem = np.divmod(gid, nx * ny)
    j, i = np.divmod(rem, nx)
    P, _ = g.cell_corners(np.stack([i, j, k], 1))
    vol = np.zeros(g.n); h = np.zeros(g.n)
    vol[gid] = np.abs(_hex_volume(P))
    h[gid] = np.maximum((P[:, :4, 2] - P[:, 4:, 2]).mean(1), 1e-6)
    P_ = md.props
    porv = np.nan_to_num(P_["PORV"]) if "PORV" in P_ else \
        np.nan_to_num(P_["PORO"]) * np.nan_to_num(P_.get("NTG", np.ones(g.n)), nan=1.0) * vol
    porv = np.where(act.ravel(), porv, 0.0)
    ntg = np.nan_to_num(P_.get("NTG", np.ones(g.n)), nan=1.0)

    def R(a):
        return np.asarray(a, float).reshape(nz, ny, nx)

    def agg(fn):
        out = []
        for k0, k1 in groups:
            out.append(fn(slice(k0, k1 + 1)))
        return np.stack(out).ravel()

    A, PV, V, H = act, R(porv), R(vol), R(h)
    newp = {}
    newp["PORV"] = agg(lambda s_: PV[s_].sum(0))
    for nm in ("PERMX", "PERMY"):
        if nm in P_:
            kk = R(np.nan_to_num(P_[nm])) * R(ntg)
            newp[nm] = agg(lambda s_: (kk[s_] * H[s_] * A[s_]).sum(0) / np.maximum((H[s_] * A[s_]).sum(0), 1e-30))
    if "PERMZ" in P_:
        kz = R(np.nan_to_num(P_["PERMZ"]))
        newp["PERMZ"] = agg(lambda s_: np.where((H[s_] * A[s_]).sum(0) > 0,
                                                (H[s_] * A[s_]).sum(0) / np.maximum(
                                                    (H[s_] * A[s_] / np.maximum(kz[s_], 1e-12)).sum(0), 1e-30), 0.0))
    pvsum = agg(lambda s_: PV[s_].sum(0))
    for nm in ("SWL", "SWATINIT", "PCW", "SWLPC", "KRW", "KRG", "IKRW", "IKRG"):
        if nm in P_:
            x = R(np.nan_to_num(P_[nm]))
            newp[nm] = agg(lambda s_: (x[s_] * PV[s_]).sum(0)) / np.maximum(pvsum, 1e-30)
    for nm in ("SATNUM", "EQLNUM", "FIPNUM", "PVTNUM", "MULTNUM", "FLUXNUM", "OPERNUM", "IMBNUM"):
        if nm in P_:
            x = R(np.nan_to_num(P_[nm])).astype(int)
            vals = np.unique(x[A])
            def maj(s_, x=x, vals=vals):
                best = np.zeros(x.shape[1:]); arg = np.full(x.shape[1:], vals[0])
                for v in vals:
                    w_ = (PV[s_] * (x[s_] == v)).sum(0)
                    arg = np.where(w_ > best, v, arg); best = np.maximum(best, w_)
                return arg
            newp[nm] = agg(maj).astype(float)
    newp["PORO"] = np.ones(nx * ny * nzc)   # porosity comes from PORV / volume
    newp["NTG"] = np.ones(nx * ny * nzc)
    actc = agg(lambda s_: (A[s_] & (PV[s_] > 0)).any(0)).astype(int)
    rst = {}
    if md.restart:
        for nm, val in md.restart.items():
            x = R(np.nan_to_num(val))
            rst[nm] = agg(lambda s_: (x[s_] * PV[s_]).sum(0)) / np.maximum(pvsum, 1e-30)
    out = Model()
    out.units, out.deck, out.mapaxes = md.units, md.deck, md.mapaxes
    out.notes = list(md.notes)
    out.grid = CornerPointGrid(nx, ny, nzc, g.coord.ravel(), zcc.ravel(), actc)
    out.props, out.restart = newp, rst
    kmap = np.zeros(nz + 1, int)
    for kc, (k0, k1) in enumerate(groups):
        kmap[k0 + 1:k1 + 2] = kc + 1
    out.kmap = kmap
    log("  vertical upscaling: %d -> %d layers (%d fine layers each), %d active cells"
        % (nz, nzc, nmerge, int(actc.sum())))
    return out


# ----------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------
def fit_corey(s, kr, s_res_other):
    """Fit MOOSE PorousFlowRelativePermeabilityCorey to a tabulated curve.
    Returns s_res, n, scaling (endpoint)."""
    s, kr = np.asarray(s, float), np.asarray(kr, float)
    o = np.argsort(s)
    s, kr = s[o], kr[o]
    zero = kr <= 1e-9
    s_res = float(s[zero].max()) if zero.any() else float(s.min())
    s_max = 1.0 - s_res_other
    scaling = float(np.interp(s_max, s, kr))
    se = (s - s_res) / max(1.0 - s_res - s_res_other, 1e-6)
    ok = (se > 1e-3) & (se < 0.999) & (kr > 1e-9) & (scaling > 0)
    if ok.sum() >= 2:
        n = float(np.sum(np.log(se[ok]) * np.log(kr[ok] / scaling)) / np.sum(np.log(se[ok]) ** 2))
    else:
        n = 2.0
    return s_res, max(n, 1.0), max(scaling, 1e-6)


def fit_brooks_corey(sw, pc):
    """Brooks-Corey fit Pc = pe * Se**(-1/lambda), Se = (Sw - x)/(1 - x), to a Pc table, with x < Swl chosen
    so the curve stays finite and smooth through the whole table, including its end point at Swl.

    Returns (swl, pc_at_swl, pe, lambda, x, rms relative error) or None if the table has no capillary pressure.
    """
    sw = np.asarray(sw, float)
    pc = np.asarray(pc, float)
    if pc.max() <= 0:
        return None
    i0 = int(np.argmin(sw))
    swl, pcmax = float(sw[i0]), float(pc[i0])
    sel = pc > 0
    if sel.sum() < 3:
        return None
    best = None
    span = 1.0 - swl
    for x in np.concatenate([swl - span * np.geomspace(1e-4, 0.5, 200)]):
        if x < 0 and swl > 0:
            x = max(x, 0.0) if swl - max(x, 0.0) > 1e-6 else None
        if x is None or x >= swl:
            continue
        se = (sw[sel] - x) / (1.0 - x)
        slope, icpt = np.polyfit(np.log(se), np.log(pc[sel]), 1)
        if slope >= 0:
            continue
        lam, pe = -1.0 / slope, float(np.exp(icpt))
        fit = pe * se ** (-1.0 / lam)
        err = float(np.sqrt(np.mean((fit / pc[sel] - 1.0) ** 2)))
        if best is None or err < best[-1]:
            best = (swl, pcmax, pe, float(lam), float(x), err)
    return best


def brooks_corey_pc(sw, pe, lam, sat_lr, pc_max):
    se = np.clip((np.asarray(sw, float) - sat_lr) / (1.0 - sat_lr), 1e-12, None)
    return np.where(se >= 1.0, pe, np.minimum(pe * se ** (-1.0 / lam), pc_max))


def _relperm_tables(relperm, m, deck, P, gi, satnum, args, rep):
    """--relperm table/deck (finite-volume input only): keep the deck's tables instead of Corey fits, and with
    SATOPTS HYSTER (--relperm deck) add each block's imbibition table (IMBNUM) and Eclipse's hysteresis option,
    for the add-on ADPorousFlowRelativePermeabilityHysteresisTable."""
    mode = getattr(args, "relperm", "corey")
    if mode == "deck" and not args.fv:
        rep("- Rel-perm: Corey fits (the tabulated / hysteresis rel-perm is only written with --fv).")
        return
    if mode == "corey" or not deck:
        return
    hy = deck.hysteresis() if mode == "deck" else None
    if mode == "deck" and hy is None:
        rep("- Rel-perm: the deck's tables (no SATOPTS HYSTER: drainage curves only).")
    if hy:
        km = hy["kr_model"]
        if km not in (2, 3):
            rep("- WARNING: EHYSTR option %d is not supported (only Killough for gas, options 2 and 3); "
                "drainage curves are used." % km)
            hy = None
        else:
            rep("- Rel-perm hysteresis (SATOPTS HYSTER, EHYSTR option %d): gas follows Killough scanning curves "
                "between the SATNUM (drainage) and IMBNUM (imbibition) tables, using each cell's largest gas "
                "saturation so far (aux variable sg_max, updated at the end of each time step); water uses the %s "
                "table%s." % (km, "imbibition" if km == 3 else "drainage",
                              " in every cell, as Eclipse does for this option" if km == 3 else ""))
            if hy["applies_to"] in ("PC", "BOTH"):
                rep("- WARNING: EHYSTR asks for capillary-pressure hysteresis (%s); not modelled, the drainage Pc "
                    "is used." % hy["applies_to"])
    imbnum = np.nan_to_num(P["IMBNUM"]).astype(int)[gi] if (hy and "IMBNUM" in P) else None
    if hy and imbnum is None:
        rep("- WARNING: SATOPTS HYSTER but no IMBNUM in the INIT file: no hysteresis.")
        hy = None
    blocks = m.blocks
    for rp in relperm:
        rel = deck.satfunc(rp["table"])
        if not rel or "sg" not in rel or "sw" not in rel:
            continue
        o_w, o_g = np.argsort(rel["sw"]), np.argsort(rel["sg"])
        tab = dict(sw=rel["sw"][o_w], krw=rel["krw"][o_w], sg=rel["sg"][o_g], krg=rel["krg"][o_g])
        bid = [b for b, nm in m.block_names.items() if nm == rp["block"]]
        sel = np.isin(blocks, bid) if bid else np.ones(len(gi), bool)

        def med(nm):
            return float(np.median(np.asarray(P[nm], float)[gi][sel])) if nm in P and sel.any() else None
        kw_cell, kg_cell = med("KRW"), med("KRG")
        tab["mult_w"] = kw_cell / tab["krw"][-1] if kw_cell and tab["krw"][-1] > 0 else 1.0
        if kg_cell and tab["krg"][-1] > 0 and abs(kg_cell / tab["krg"][-1] - 1) > 0.01:
            rep("- NOTE: %s: KRG %.4g differs from the table end point %.4g (not scaled)."
                % (rp["block"], kg_cell, tab["krg"][-1]))
        if hy:
            vals, cnts = np.unique(imbnum[sel], return_counts=True)
            it = int(vals[np.argmax(cnts)])
            ri = deck.satfunc(it)
            if ri and "sg" in ri and "sw" in ri:
                oi_w, oi_g = np.argsort(ri["sw"]), np.argsort(ri["sg"])
                tab.update(imb=it, isw=ri["sw"][oi_w], ikrw=ri["krw"][oi_w], isg=ri["sg"][oi_g], ikrg=ri["krg"][oi_g],
                           kr_model=hy["kr_model"])
                ik = med("IKRW")
                tab["mult_wi"] = ik / tab["ikrw"][-1] if ik and tab["ikrw"][-1] > 0 else 1.0
                if len(vals) > 1:
                    rep("- NOTE: %s has cells of IMBNUM %s; table %d (the most common) is used."
                        % (rp["block"], " ".join(map(str, vals.tolist())), it))
        rp["tab"] = tab
    done = set()
    for rp in relperm:
        t = rp.get("tab")
        if not t or rp["table"] in done:
            continue
        done.add(rp["table"])
        line = "- Rel-perm table %d (SWFN/SGFN, linear interpolation)" % rp["table"]
        if "imb" in t:
            line += ", imbibition table %d: trapped-gas (critical) saturation %.4f (drainage %.4f)" % (
                t["imb"], float(t["isg"][t["ikrg"] <= 0].max()) if (t["ikrg"] <= 0).any() else t["isg"][0],
                float(t["sg"][t["krg"] <= 0].max()) if (t["krg"] <= 0).any() else t["sg"][0])
            if t["kr_model"] == 3:
                line += "; water krw at the end of the imbibition table %.4g scaled by IKRW to %.4g (x%.3g)" % (
                    t["ikrw"][-1], t["ikrw"][-1] * t["mult_wi"], t["mult_wi"])
        rep(line + ".")


def conn_wells(p, args):
    """Wells written with the connection model (--fv --well-model connection, gas production with COMPDAT CF)."""
    if not args.fv or getattr(args, "well_model", "rate") != "connection":
        return []
    return [wl for wl in p["wells"] if wl.get("cf_var") and any(q < 0 for q in wl["qg"])]


def _lst(a, nd=5):
    return "'%s'" % " ".join(("%%.%dg" % nd) % float(v) for v in a)


def _write_relperm_tables(blk, rp, multi):
    """Tabulated rel-perm of both phases for one block (add-on ADPorousFlowRelativePermeabilityHysteresisTable)."""
    t = rp["tab"]
    hy = "imb" in t
    # gas: Eclipse two-point end-point scaling of the Sg axis to SGU = 1 - SWL of the block's cells (as the
    # Corey fits do with swr_gas); the imbibition table keeps its critical saturation
    f = 1.0
    if rp.get("swr_gas") is not None and t["sg"][-1] > 0:
        f = (1.0 - rp["swr_gas"]) / t["sg"][-1]
    sg = t["sg"] * f
    for ph, nm in ((0, "water"), (1, "gas")):
        prm = [("block", rp["block"])] if multi else []
        prm.append(("phase", ph))
        if ph == 0:
            prm += [("saturation", _lst(t["sw"])), ("kr", _lst(t["krw"]))]
            if hy and t["kr_model"] == 3:
                prm += [("imbibition_saturation", _lst(t["isw"])), ("imbibition_kr", _lst(t["ikrw"])),
                        ("hysteresis", "imbibition # EHYSTR option 3: imbibition curve for water in every cell"),
                        ("kr_multiplier", "%.6g # IKRW / last kr of IMBNUM table %d" % (t["mult_wi"], t["imb"]))]
                cm = "water: SATNUM table %d (drainage), IMBNUM table %d (imbibition, used)" % (rp["table"], t["imb"])
            else:
                if abs(t["mult_w"] - 1) > 1e-6:
                    prm.append(("kr_multiplier", "%.6g # KRW / last kr of the table" % t["mult_w"]))
                cm = "water: SATNUM table %d" % rp["table"]
        else:
            prm += [("saturation", _lst(sg)), ("kr", _lst(t["krg"]))]
            cm = "gas: SATNUM table %d%s" % (rp["table"], "" if f == 1.0 else ", Sg axis x %.4f (SGU = 1 - SWL)" % f)
            if hy:
                sc = float(t["isg"][t["ikrg"] <= 0].max()) if (t["ikrg"] <= 0).any() else float(t["isg"][0])
                g = (f * t["isg"][-1] - sc) / max(t["isg"][-1] - sc, 1e-12)
                isg = sc + (t["isg"] - sc) * g
                prm += [("imbibition_saturation", _lst(isg)), ("imbibition_kr", _lst(t["ikrg"])),
                        ("hysteresis", "killough"), ("max_saturation", "sg_max")]
                cm += "; Killough scanning curves towards IMBNUM table %d" % t["imb"]
        blk("relperm_%s%s" % (nm, "_" + rp["block"] if multi else ""), "ADPorousFlowRelativePermeabilityHysteresisTable",
            prm, cm)


def _capillary_blocks(m, relperm, deck, P, gi, sg_e, restart, args):
    """Split each rel-perm block into sub-blocks with their own Brooks-Corey capillary pressure.

    Eclipse scales each cell's Pc table to the cell's PCW (Pc at SWLPC, honouring SWATINIT) and
    SWLPC. MOOSE takes one Pc curve per block, so cells are binned on log(PCW) and SWLPC.
    """
    lines = []
    fits = {}
    for rp in relperm:
        t = rp["table"]
        if t in fits:
            continue
        rel = deck.satfunc(t) if deck else None
        fits[t] = fit_brooks_corey(rel["sw"], rel["pc"]) if rel and "pc" in rel else None
        f = fits[t]
        if f:
            lines.append("- Capillary pressure table %d: Brooks-Corey fit pe=%.4g Pa lambda=%.3f, Se=(Sw-%.4f)/(1-%.4f) "
                         "(table Swl=%.3f, Pc(Swl)=%.4g Pa; rms error %.1f%% over the table)."
                         % (t, f[2], f[3], f[4], f[4], f[0], f[1], 100 * f[5]))
        else:
            lines.append("- Capillary pressure table %d: zero, kept at Pc = 0." % t)
    n = len(gi)
    swl_c = np.asarray(P["SWL"], float)[gi] if "SWL" in P else np.full(n, np.nan)
    pcw = np.nan_to_num(P["PCW"])[gi] if "PCW" in P else np.full(n, np.nan)
    swlpc = np.nan_to_num(P["SWLPC"], nan=np.nan)[gi] if "SWLPC" in P else np.full(n, np.nan)
    by_name = {rp["block"]: rp for rp in relperm}
    new_blocks = m.blocks.copy()
    new_names = {}
    new_rel = []
    nid = 0
    for b_ in sorted(set(m.blocks.tolist())):
        base = by_name[m.block_names[b_]]
        f = fits.get(base["table"])
        ids = np.nonzero(m.blocks == b_)[0]
        if f is None:
            new_blocks[ids] = nid
            new_names[nid] = base["block"]
            new_rel.append(dict(base))
            nid += 1
            continue
        swl_t, pcmax_t, pe_t, lam, x_t, _ = f
        c_pcw = np.where(np.isfinite(pcw[ids]) & (pcw[ids] > 0), pcw[ids], pcmax_t)
        c_swl = np.where(np.isfinite(swlpc[ids]), swlpc[ids], swl_t)
        kp = np.round(np.log(c_pcw) / np.log(args.pc_bin_ratio)).astype(int)
        ks = np.round(c_swl / args.pc_swl_bin).astype(int)
        keys = sorted(set(zip(kp.tolist(), ks.tolist())))
        for q, (a, b) in enumerate(keys):
            sel = ids[(kp == a) & (ks == b)]
            sub = (kp == a) & (ks == b)
            pcw_rep = float(np.exp(np.mean(np.log(c_pcw[sub]))))
            swl_rep = float(np.clip(np.mean(c_swl[sub]), 0.0, 0.95))
            # Eclipse scaling of the table's Sw axis [Swl_t, 1] onto [SWLPC, 1] keeps the Brooks-Corey form,
            # with the shifted residual mapped the same way
            sat_lr = float(np.clip(swl_rep - (swl_t - x_t) * (1.0 - swl_rep) / (1.0 - swl_t), 0.0, 0.95))
            nm = base["block"] if len(keys) == 1 else "%s_pc%02d" % (base["block"], q + 1)
            new_blocks[sel] = nid
            new_names[nid] = nm
            r2 = dict(base)
            # gas rel-perm end point at the cells' SWL (Eclipse end-point scaling: SGU = 1 - SWL), so the gas
            # curve has no kink inside the saturations the gas cap actually reaches
            swl_g = swl_c[ids][sub]
            if np.isfinite(swl_g).any():
                r2["swr_gas"] = float(np.clip(np.nanmean(swl_g), 0.0, 0.95))
            r2.update(block=nm, pc=dict(pe=pe_t * pcw_rep / pcmax_t, lam=lam, sat_lr=sat_lr,
                                        pc_max=args.pc_max_factor * pcw_rep, pcw=pcw_rep))
            new_rel.append(r2)
            nid += 1
    m.blocks = new_blocks
    m.block_names = new_names
    # equilibrium check: model Pc at the initial Sw against Eclipse's PCOW
    if sg_e is not None and restart and "PCOW" in restart:
        pc_mod = np.zeros(n)
        for b_, nm in new_names.items():
            rp = [r for r in new_rel if r["block"] == nm][0]
            sel = m.blocks == b_
            if "pc" in rp:
                c = rp["pc"]
                pc_mod[sel] = brooks_corey_pc(1.0 - sg_e[sel], c["pe"], c["lam"], c["sat_lr"], c["pc_max"])
        pc_ecl = restart["PCOW"][gi]
        gas = sg_e > 1e-3
        if gas.any():
            d = (pc_mod - pc_ecl)[gas] / 1e5
            lines.append("- Capillary pressure: %d blocks (cells binned by PCW ratio %.2f and SWLPC %.3f). At the "
                         "initial saturations of the %d gas-bearing cells, MOOSE Pc - Eclipse PCOW: median %.3f bar, "
                         "95th pct |.| %.3f bar (Eclipse PCOW there: median %.3f bar). The immobile-water surrogate "
                         "is off: water above residual is held by capillarity."
                         % (len(new_names), args.pc_bin_ratio, args.pc_swl_bin, int(gas.sum()), np.median(d),
                            np.percentile(np.abs(d), 95), np.median(pc_ecl[gas]) / 1e5))
    return new_rel, lines


def quad_area(p):
    """p: (m,4,3) quad corners in cyclic order -> area"""
    d1, d2 = p[:, 2] - p[:, 0], p[:, 3] - p[:, 1]
    return 0.5 * np.linalg.norm(np.cross(d1, d2), axis=1)


_SIDE_NODES = {1: [0, 1, 5, 4], 2: [1, 2, 6, 5], 3: [2, 3, 7, 6], 4: [3, 0, 4, 7], 5: [0, 3, 2, 1], 6: [4, 5, 6, 7]}


def sideset_area(m, name):
    e, s = m.sidesets[name]
    if len(e) == 0:
        return 0.0
    idx = np.array([_SIDE_NODES[int(x)] for x in s])
    pts = m.nodes[np.take_along_axis(m.conn[e], idx, 1)]
    return float(quad_area(pts).sum())


def safe(name):
    out = "".join(c if c.isalnum() or c == "_" else "_" for c in str(name))
    return out if out and not out[0].isdigit() else "w_" + out


def fmt(x):
    return "%.10g" % (x + 0.0)



def rate_sync_times(wells, end):
    """Times the time stepper must hit: the schedule knots next to a rate change (both sides), so every
    production or injection spell is resolved, but not the report dates where nothing changes (after a
    shut-in Eclipse reports daily, which would otherwise force one-day steps for years)."""
    keep = set()
    for wl in wells:
        t = list(wl["t"])
        for q in (wl.get("qg", []), wl.get("qw", [])):
            q = list(q)
            for i in range(len(t)):
                if (i > 0 and q[i] != q[i - 1]) or (i + 1 < len(t) and q[i] != q[i + 1]):
                    keep.add(t[i])
    return sorted(x for x in keep if 0 < x < end)

def wrap_list(vals, per=8, indent="        "):
    vals = [fmt(v) for v in vals]
    rows = [" ".join(vals[i:i + per]) for i in range(0, len(vals), per)]
    return ("\n" + indent).join(rows)


# ----------------------------------------------------------------------------------------------
# main conversion
# ----------------------------------------------------------------------------------------------
def convert(args, log=print):
    report = []
    rep = report.append

    # --- load geometry + properties
    deck = None
    if args.case:
        deck = eclio.Deck(args.case)
        log("  DATA deck %s: units %s, phases %s" % (os.path.basename(args.case), deck.units,
                                                     ",".join(sorted(deck.phases)) or "?"))
        args.egrid = args.egrid or _find_case_file(args.case, ".EGRID")
        args.init = args.init or _find_case_file(args.case, ".INIT")
        if args.restart_step is not None:
            args.unrst = args.unrst or _find_case_file(args.case, ".UNRST")
            if not args.unrst:  # non-unified restart files CASE.Xnnnn
                base = os.path.splitext(args.case)[0]
                d = os.path.dirname(os.path.abspath(args.case))
                stem = os.path.basename(base).lower()
                xs = sorted(f for f in os.listdir(d) if f.lower().startswith(stem + ".x") and f[-4:].isdigit())
                if xs:
                    pick = xs[-1] if args.restart_step < 0 else next(
                        (f for f in xs if int(f[-4:]) == args.restart_step), None)
                    if pick:
                        args.unrst = os.path.join(d, pick)
                        args.restart_step = 0
        pet = deck.get("PETGRID")
        if not args.egrid and not args.grdecl and pet and pet[0]:
            cand = os.path.join(os.path.dirname(os.path.abspath(args.case)), str(pet[0][0]))
            try:
                args.gsg = eclio._resolve_path(os.path.dirname(os.path.abspath(args.case)), str(pet[0][0]))
            except FileNotFoundError:
                raise SystemExit("Deck uses PETGRID %s but the file is missing" % cand)
            rep("- Grid geometry decoded from Petrel's PETGRID file %s (no EGRID was written: GRIDFILE 0 0)."
                % os.path.basename(args.gsg))
        if not args.egrid and not args.grdecl and not getattr(args, "gsg", None):
            # case exported but never run: the grid is in the deck's GRID-section (GRDECL) includes
            log("  no .EGRID next to the deck: reading COORD/ZCORN/ACTNUM/properties from the deck itself")
            args.grdecl = args.case
            report.append("- Grid and properties read from the deck's GRDECL includes (no EGRID/INIT). "
                          "EQUALS/COPY/MULTIPLY/BOX edits in the deck are NOT applied - check porosity/permeability.")
    if getattr(args, "gsg", None):
        md = load_from_gsg(args.gsg, args.init, args.unrst if args.restart_step is not None else None,
                           args.restart_step or 0, log)
        for note in md.notes:
            rep("- " + note)
    elif args.egrid:
        md = load_from_egrid(args.egrid, args.init, args.unrst if args.restart_step is not None else None,
                             args.restart_step or 0, log)
    elif args.grdecl:
        md = load_from_grdecl(args.grdecl, args.props or (), log)
    else:
        raise SystemExit("Need --case, --egrid or --grdecl")
    if args.props and args.egrid:
        for p in args.props:
            for k, a in eclio.read_grdecl(p).items():
                if len(a) == md.grid.n and a.dtype.kind == "f":
                    md.props[k] = a
    if args.egrid and args.case and not args.init:
        # EGRID but no INIT: take properties from the deck's GRDECL includes
        for k, a in eclio.read_grdecl(args.case).items():
            if len(a) == md.grid.n and a.dtype.kind == "f" and k not in ("ZCORN", "ACTNUM") and k not in md.props:
                md.props[k] = a
        report.append("- No .INIT found: properties read from the deck's GRDECL includes (deck edits not applied).")
    md.deck = deck
    if deck and not args.init:
        md.units = deck.units
    md.kmap = None
    if args.merge_layers and args.merge_layers > 1:
        md = upscale_layers(md, args.merge_layers, log)
        rep("- Vertical upscaling: every %d layers merged (%d -> %d). Pore volume and restart saturations are "
            "PV-conserved, kh arithmetic, kv harmonic, regions by PV majority; wells and aquifers remapped."
            % (args.merge_layers, int(md.kmap.max() and len(md.kmap) - 1), md.grid.nz))
    grid, P = md.grid, md.props
    n = grid.n

    # --- properties (SI)
    if "PORO" not in P:
        raise SystemExit("No PORO found (need .INIT, or PORO in GRDECL/--props)")
    if "PERMX" not in P:
        raise SystemExit("No PERMX found (need .INIT, or PERMX in GRDECL/--props)")
    poro = np.nan_to_num(P["PORO"])
    ntg = np.nan_to_num(P.get("NTG", np.ones(n)), nan=1.0)
    multpv = np.nan_to_num(P.get("MULTPV", np.ones(n)), nan=1.0)
    phi = poro * ntg * multpv
    if "NTG" in P and "PORV" not in P:
        rep("- Porosity passed to MOOSE is PORO x NTG%s (Eclipse pore-volume convention)."
            % (" x MULTPV" if "MULTPV" in P else ""))
    kx = np.nan_to_num(P["PERMX"])
    ky = np.nan_to_num(P.get("PERMY", kx))
    if "PERMZ" in P:
        kz = np.nan_to_num(P["PERMZ"])
    else:
        kz = kx * args.kvkh
        rep("- PERMZ missing: used PERMX x %g." % args.kvkh)
    if "PERMY" not in P:
        rep("- PERMY missing: used PERMX.")
    # Eclipse PERMX/PERMY are net-to-gross corrected in transmissibility for X/Y only
    if "NTG" in P and not args.no_ntg_perm:
        kx, ky = kx * ntg, ky * ntg
        rep("- Horizontal permeability multiplied by NTG (Eclipse applies NTG to X/Y transmissibility).")

    phi_from_porv = "PORV" in P and getattr(args, "use_porv", True)
    keep = (phi > args.min_poro) | (phi_from_porv & (np.nan_to_num(P.get("PORV", np.zeros(n))) > 0))
    keep &= phi > 0
    rep("- Cells with effective porosity <= %g are excluded: %d active cells." %
        (args.min_poro, int(((grid.actnum.ravel() > 0) & ~keep).sum())))

    # --- mesh
    log("Building mesh ...")
    split = None
    if deck:
        seals = [(r1, r2, a) for r1, r2, mult, a in deck.multregt() if mult == 0.0]
        if seals:
            arr = {"M": "MULTNUM", "F": "FLUXNUM", "O": "OPERNUM"}.get(seals[0][2], "MULTNUM")
            if arr in P:
                split = np.nan_to_num(P[arr]).astype(int)
                rep("- MULTREGT zero-transmissibility barrier(s) between %s regions %s: the mesh is split so "
                    "no nodes are shared across them (sealing)." % (arr, ", ".join("%d-%d" % (a, b) for a, b, _ in seals)))
            else:
                rep("- WARNING: MULTREGT barrier uses %s, which is not in the INIT/GRDECL; barrier ignored." % arr)
    use_nnc = bool(args.fv and args.nnc == "on" and args.egrid and args.init
                   and os.path.exists(args.egrid) and os.path.exists(args.init))
    if args.fault_snap is None:
        args.fault_snap = 0.0 if use_nnc else 1.0
    m = build_mesh(grid, keep=keep, merge_tol=args.merge_tol, shift_origin=not args.no_shift,
                   split_region=split, log=log, fault_snap=args.fault_snap, fault_snap_frac=args.fault_snap_frac)
    gi = m.gidx
    pv_all = (phi * 1.0)[(grid.actnum.ravel() > 0) & keep]
    m.info = ["petrel2moose conversion", "origin_x %.3f origin_y %.3f (add to x,y for map coords)" % tuple(m.origin[:2]),
              "z is elevation in m (negative depth)", "perm in m2, pressure in Pa"]
    rep("- Mesh: %d HEX8 elements, %d nodes, %d elements with mirrored node order (left-handed grid)."
        % (m.stats["n_elem"], m.stats["n_nodes"], m.stats["n_flipped"]))
    if args.fault_snap > 0:
        rep("- Faults with small throw: %d nodes merged across faults whose throw is below min(%g m, %g x the cell "
            "thickness), so those faces connect as in Eclipse (they would otherwise be no-flow)."
            % (m.stats.get("fault_snapped_nodes", 0), args.fault_snap, args.fault_snap_frac))
    if not use_nnc:
        rep("- %d of %d logical internal faces are non-conforming (fault throw); MOOSE treats them as no-flow."
            % (m.stats["nonconforming_faces"], m.stats["internal_faces"]))
    if md.mapaxes is not None:
        rep("- MAPAXES present (%s); coordinates are used as stored in COORD (no transform applied)."
            % " ".join("%.1f" % v for v in md.mapaxes))
    rep("- Coordinates shifted by origin (%.1f, %.1f) m; z = -depth." % tuple(m.origin[:2]))

    # --- permeability tensor rotated into global xy (PERMX is along the grid I axis)
    nodes_e = m.nodes[m.conn]          # (E,8,3), Exodus order
    loc = np.where(m.flipped[:, None], np.array([0, 3, 2, 1, 4, 7, 6, 5])[None, :], np.arange(8)[None, :])
    nl = np.take_along_axis(nodes_e, loc[:, :, None], 1)  # un-flipped local order
    eI = nl[:, [1, 2, 5, 6]].mean(1) - nl[:, [0, 3, 4, 7]].mean(1)
    ang = np.arctan2(eI[:, 1], eI[:, 0])
    c, s = np.cos(ang), np.sin(ang)
    KX, KY, KZ = kx[gi] * eclio.MD_TO_M2, ky[gi] * eclio.MD_TO_M2, kz[gi] * eclio.MD_TO_M2
    m.elem_vars["poro"] = phi[gi]
    m.elem_vars["perm_xx"] = c * c * KX + s * s * KY
    m.elem_vars["perm_yy"] = s * s * KX + c * c * KY
    m.elem_vars["perm_xy"] = c * s * (KX - KY)
    m.elem_vars["perm_zz"] = KZ
    m.elem_vars["permx_md"] = kx[gi]
    m.elem_vars["ntg"] = ntg[gi]
    m.elem_vars["cell_i"] = m.cells[:, 0] + 1.0
    m.elem_vars["cell_j"] = m.cells[:, 1] + 1.0
    m.elem_vars["cell_k"] = m.cells[:, 2] + 1.0
    for reg in ("SATNUM", "PVTNUM", "FIPNUM", "EQLNUM", "MULTNUM"):
        if reg in P:
            m.elem_vars[reg.lower()] = np.nan_to_num(P[reg])[gi]
    vol = np.abs(_hex_volume(nodes_e))
    if phi_from_porv:
        pv_si = np.nan_to_num(P["PORV"])[gi] * eclio.UNITS[md.units]["liq_rate"]  # rm3/rb -> m3
        phi_porv = np.clip(pv_si / vol, 0.0, 1.0)
        if getattr(md, "kmap", None) is not None:
            # upscaled model: PORO was replaced by 1 (porosity comes from PORV), so no ratio to report
            rep("- Porosity set from Eclipse PORV / cell volume so pore volume matches Eclipse exactly "
                "(upscaled cells: median porosity %.4f, 1st-99th pct %.3f-%.3f)."
                % (np.median(phi_porv), np.percentile(phi_porv, 1), np.percentile(phi_porv, 99)))
        else:
            ratio = phi_porv / np.maximum(phi[gi], 1e-12)
            rep("- Porosity set from Eclipse PORV / cell volume so pore volume matches Eclipse exactly "
                "(median ratio to PORO*NTG*MULTPV: %.4f, 1st-99th pct %.3f-%.3f)."
                % (np.median(ratio), np.percentile(ratio, 1), np.percentile(ratio, 99)))
        m.elem_vars["poro"] = phi_porv
    m.elem_vars["pore_volume"] = vol * m.elem_vars["poro"]
    if phi_from_porv:
        pv_ecl = float(np.nansum(np.nan_to_num(P["PORV"])[grid.actnum.ravel() > 0])) * eclio.UNITS[md.units]["liq_rate"]
    else:
        pv_ecl = _eclipse_pv(md, phi, keep)
    rep("- Pore volume in MOOSE mesh: %.6g m3 (%.4f%% of Eclipse active pore volume kept)."
        % (m.elem_vars["pore_volume"].sum(), 100.0 * m.elem_vars["pore_volume"].sum() / max(pv_ecl, 1e-30)))

    # --- blocks
    m.blocks = np.zeros(len(gi), int)
    m.block_names = {0: "reservoir"}
    if args.zones:
        for bid, z in enumerate(args.zones.split(","), 1):
            name, rng = z.split(":")
            k1, k2 = [int(v) for v in rng.split("-")]
            m.blocks[(m.cells[:, 2] + 1 >= k1) & (m.cells[:, 2] + 1 <= k2)] = bid
            m.block_names[bid] = safe(name)
    elif args.block_by:
        key = args.block_by.upper()
        if key not in P:
            raise SystemExit("--block-by %s: property not found" % key)
        vals = np.nan_to_num(P[key])[gi].astype(int)
        m.blocks = vals
        m.block_names = {v: "%s_%d" % (key.lower(), v) for v in set(vals.tolist())}
    used = set(m.blocks.tolist())
    m.block_names = {b: nm for b, nm in m.block_names.items() if b in used}

    # --- fluid & rock data from deck
    rho_o, rho_w, rho_g_sc = (deck.densities() if deck else (859.0, 1030.0, 0.68))
    rho_w_sc = rho_w
    pvtw = deck.pvtw() if deck else None
    pvdg = deck.pvdg() if deck else None
    if pvtw:
        rho_w = rho_w_sc / pvtw["bw"]        # reservoir-condition water density for hydrostatics
        rep("- Water PVT from PVTW: Bw=%.4f at %.4g Pa, cw=%.3g 1/Pa, mu=%.3g Pa s -> SimpleFluidProperties "
            "(reservoir density %.1f kg/m3)." % (pvtw["bw"], pvtw["pref"], pvtw["cw"], pvtw["mu"], rho_w))
    gas_pvt_file = None
    if pvdg:
        pg, bg, mug = pvdg
        o_ = np.argsort(pg)
        pg, bg, mug = pg[o_], bg[o_], mug[o_]
        inv_bg = 1.0 / bg

        def rho_g_of_p(p):
            p = np.asarray(p, float)
            # linear in 1/Bg; through the origin below the table, linear extrapolation above it
            x = np.concatenate([[0.0], pg])
            y = np.concatenate([[0.0], inv_bg])
            slope = (y[-1] - y[-2]) / (x[-1] - x[-2])
            val = np.where(p <= x[-1], np.interp(p, x, y), y[-1] + slope * (p - x[-1]))
            return rho_g_sc * val

        def mu_g_of_p(p):
            p = np.asarray(p, float)
            slope = (mug[-1] - mug[-2]) / (pg[-1] - pg[-2])
            return np.where(p <= pg[-1], np.interp(p, pg, mug), mug[-1] + slope * (p - pg[-1]))

        gas_pvt_file = args.output + "_gas_pvt.csv"
        pgrid = np.unique(np.concatenate([np.linspace(1e5, max(6e7, 1.5 * pg.max()), 300), pg]))
        tgrid = np.array([273.15, 323.15, 373.15, 423.15, 473.15])
        with open(gas_pvt_file, "w") as f:
            f.write("pressure,temperature,density,viscosity\n")
            for p_ in pgrid:
                r_, m_ = float(rho_g_of_p(p_)), float(mu_g_of_p(p_))
                for t_ in tgrid:
                    f.write("%.8g,%.8g,%.8g,%.8g\n" % (p_, t_, r_, m_))
        rep("- Gas PVT from PVDG (%d points, %.3g-%.3g Pa) tabulated to %s for TabulatedBicubicFluidProperties "
            "(isothermal: same values at every T)." % (len(pg), pg.min(), pg.max(), os.path.basename(gas_pvt_file)))
    else:
        rho_g_of_p = None
    T = (deck.temperature() if deck else None) or args.temperature + 273.15
    if (not deck or deck.temperature() is None) and not (deck and deck.pvdg() and deck.pvtw()):
        rep("- Reservoir temperature not in deck (RTEMP): using %.1f C (--temperature)." % (T - 273.15))

    # relative permeability: one Corey set per block; blocks follow SATNUM when there are several tables
    ntab = deck.n_satfunc() if deck else 0
    satnum = np.nan_to_num(P["SATNUM"]).astype(int) if "SATNUM" in P else np.ones(n, int)
    sat_used = sorted(set(satnum[gi].tolist()))
    if ntab > 1 and len(sat_used) > 1 and not args.zones and not args.block_by:
        m.blocks = satnum[gi].copy()
        m.block_names = {v: "satnum_%d" % v for v in sat_used}
        rep("- %d SATNUM regions in use -> one MOOSE block per SATNUM (%s), each with its own Corey fit."
            % (len(sat_used), ", ".join(m.block_names.values())))
    relperm = []   # list of dict(block_names, swr, sgr, nw, ng, kw_max, kg_max)
    blk_ids = sorted(set(m.blocks.tolist()))
    for b in blk_ids:
        tab_no = 1
        if m.block_names.get(b, "").startswith("satnum_"):
            tab_no = b
        elif "SATNUM" in P:
            vals, cnts = np.unique(satnum[gi][m.blocks == b], return_counts=True)
            tab_no = int(vals[np.argmax(cnts)])
        rel = deck.satfunc(tab_no) if deck else None
        d_ = dict(block=m.block_names.get(b, "block_%d" % b), table=tab_no, swr=0.2, sgr=0.05, nw=3.0, ng=2.0,
                  kw_max=1.0, kg_max=1.0)
        if rel and "sw" in rel and "sg" in rel:
            swr0 = float(rel["sw"][rel["krw"] <= 1e-9].max()) if (rel["krw"] <= 1e-9).any() else float(rel["sw"].min())
            sgr0 = float(rel["sg"][rel["krg"] <= 1e-9].max()) if (rel["krg"] <= 1e-9).any() else 0.0
            d_["swr"], d_["nw"], d_["kw_max"] = fit_corey(rel["sw"], rel["krw"], sgr0)
            d_["sgr"], d_["ng"], d_["kg_max"] = fit_corey(rel["sg"], rel["krg"], swr0)
            rep("- Rel-perm %s <- Corey fit to %s: Swr=%.3f nw=%.2f krw_max=%.3f, Sgr=%.3f ng=%.2f krg_max=%.3f."
                % (d_["block"], rel["src"], d_["swr"], d_["nw"], d_["kw_max"], d_["sgr"], d_["ng"], d_["kg_max"]))
        else:
            rep("- No gas/water saturation table for %s: default Corey curves - EDIT the input." % d_["block"])
        relperm.append(d_)
    _relperm_tables(relperm, m, deck, P, gi, satnum, args, rep)
    if len(blk_ids) > 1 and len(relperm) and "satnum" not in relperm[0]["block"] and len(sat_used) > 1:
        rep("- WARNING: blocks do not follow SATNUM; each block uses its majority SATNUM table.")
    if "SWL" in P and np.nanstd(P["SWL"][gi]) > 1e-3:
        rep("- NOTE: the deck uses end-point scaling (SWL varies %.2f-%.2f cell by cell); MOOSE Corey curves use one "
            "set of end points per block. Initial saturations still honour the cell values (from the restart)."
            % (np.nanmin(P["SWL"][gi]), np.nanmax(P["SWL"][gi])))
    swr = float(np.median([r_["swr"] for r_ in relperm]))

    # --- initial state (always written as nodal fields in the mesh -> initial_from_file_var)
    cen = cell_centres(m)
    depth_e = -cen[:, 2]
    init_mode = "none"
    p_e = sg_e = None
    eqs = deck.equil_all() if deck else []
    rho_g_res = None
    for e_ in eqs:
        e_["rho_g"] = (float(rho_g_of_p(e_["datum_p"])) if rho_g_of_p else
                       e_["datum_p"] * M_CH4 / (args.zfactor * R_GAS * T))
        gwc = e_["contact1"] if e_["contact1"] is not None else e_["datum_depth"]
        if e_["datum_depth"] <= gwc:
            p_gwc = e_["datum_p"] + e_["rho_g"] * G * (gwc - e_["datum_depth"])
        else:
            p_gwc = e_["datum_p"] - rho_w * G * (e_["datum_depth"] - gwc)
        e_.update(gwc=gwc, p_gwc=p_gwc)
    eq = eqs[0] if eqs else None
    if eq:
        rho_g_res = eq["rho_g"]
    eqlnum = np.nan_to_num(P.get("EQLNUM", np.ones(n))).astype(int)

    def equil_eval(depth, region, swl):
        """Hydrostatic p and gas saturation at depths for EQUIL regions (1-based)."""
        p = np.zeros_like(depth)
        sg = np.zeros_like(depth)
        for r_, e_ in enumerate(eqs, 1):
            sel = region == r_
            if not sel.any():
                continue
            dd = depth[sel] - e_["gwc"]
            p[sel] = e_["p_gwc"] + np.where(dd < 0, e_["rho_g"], rho_w) * G * dd
            sg[sel] = np.where(dd < 0, 1.0 - swl[sel], 0.0)
        return p, sg

    swl_cell = np.nan_to_num(P["SWL"], nan=swr) if "SWL" in P else np.full(n, swr)
    def equil_nodes(reg_e):
        """EQUIL pressure/saturation evaluated at every element corner, averaged at nodes."""
        num_p = np.zeros(len(m.nodes)); num_s = np.zeros(len(m.nodes)); cnt = np.zeros(len(m.nodes))
        for c_ in range(8):
            nid = m.conn[:, c_]
            pc, sc = equil_eval(-m.nodes[nid, 2], reg_e, swl_cell[gi])
            np.add.at(num_p, nid, pc); np.add.at(num_s, nid, sc); np.add.at(cnt, nid, 1)
        return num_p / cnt, num_s / cnt

    if md.restart.get("PRESSURE") is not None:
        sg_e = md.restart["SGAS"][gi] if "SGAS" in md.restart else 1.0 - md.restart.get("SWAT", np.ones(n))[gi]
        w = m.elem_vars["pore_volume"]
        sg_n = element_to_node(m, sg_e, w)
        init_mode = "restart"
        if eqs and not args.restart_pressure:
            # Eclipse's PRESSURE is the (gas-)phase pressure including cell-scaled capillary pressure;
            # with Pc = 0 in MOOSE it is not in equilibrium. Use hydrostatic EQUIL pressure instead
            # and keep the restart saturations (so gas in place is preserved).
            reg_e = np.clip(eqlnum[gi], 1, len(eqs))
            p_e, _ = equil_eval(depth_e, reg_e, swl_cell[gi])
            p_n, _ = equil_nodes(reg_e)
            # nodal saturation must respect the contact: no gas below the GWC, and above it average only
            # cells that are themselves above the contact (avoids smearing water/gas across the GWC)
            gwc_e = np.array([eqs[r_ - 1]["gwc"] for r_ in reg_e])
            above_e = depth_e < gwc_e
            num = np.zeros(len(m.nodes)); den = np.zeros(len(m.nodes)); gwc_n = np.zeros(len(m.nodes)); cnt = np.zeros(len(m.nodes))
            for c_ in range(8):
                nid = m.conn[:, c_]
                np.add.at(num, nid, np.where(above_e, sg_e * w, 0.0))
                np.add.at(den, nid, np.where(above_e, w, 0.0))
                np.add.at(gwc_n, nid, gwc_e); np.add.at(cnt, nid, 1)
            gwc_n /= cnt
            depth_n = -m.nodes[:, 2]
            sg_n = np.where((depth_n < gwc_n) & (den > 0), num / np.maximum(den, 1e-30), 0.0)
            if args.capillary == "deck" and "PCOW" in md.restart:
                # with Eclipse's capillary pressure in MOOSE, the restart is the equilibrium state: Eclipse's
                # PRESSURE is the gas pressure, so the water pressure is PRESSURE - PCOW
                p_e = md.restart["PRESSURE"][gi] - md.restart["PCOW"][gi]
                p_n = element_to_node(m, p_e, w)
                dp = (p_e - equil_eval(depth_e, reg_e, swl_cell[gi])[0]) / 1e5
                rep("- Initial saturations from restart %s; initial water pressure = restart PRESSURE - PCOW (Eclipse "
                    "reports the gas pressure). Difference to hydrostatic EQUIL water pressure: median %.2f bar, "
                    "99th pct |.| %.2f bar." % (os.path.basename(args.unrst or ""), np.median(dp),
                                                np.percentile(np.abs(dp), 99)))
            else:
                dp = (p_e - md.restart["PRESSURE"][gi]) / 1e5
                rep("- Initial saturations from restart %s; initial pressure hydrostatic from EQUIL (MOOSE has Pc = 0, "
                    "so the restart's capillary-laden phase pressure would not be in equilibrium). Difference to the "
                    "restart pressure: median %.2f bar, 99th pct |.| %.2f bar. Use --restart-pressure to take it as is."
                    % (os.path.basename(args.unrst or ""), np.median(dp), np.percentile(np.abs(dp), 99)))
        else:
            p_e = md.restart["PRESSURE"][gi]
            p_n = element_to_node(m, p_e, w)
            rep("- Initial pressure/saturation from restart %s (pore-volume averaged to nodes)."
                % os.path.basename(args.unrst or ""))
    elif eqs:
        reg_e = np.clip(eqlnum[gi], 1, len(eqs))
        p_e, sg_e = equil_eval(depth_e, reg_e, swl_cell[gi])
        p_n, sg_n = equil_nodes(reg_e)
        init_mode = "restart"   # nodal fields in the mesh
        for r_, e_ in enumerate(eqs, 1):
            rep("- EQUIL region %d: datum %.1f m, %.4g Pa; GWC %.1f m; gas density ~%.0f kg/m3 (Z=%.2f)%s."
                % (r_, e_["datum_depth"], e_["datum_p"], e_["gwc"], e_["rho_g"], args.zfactor,
                   "" if len(eqs) == 1 else ", %d elements" % int((reg_e == r_).sum())))
        if "SWATINIT" in P:
            rep("  NOTE: the deck uses SWATINIT (capillary-scaled initial water). For an exact initial state "
                "use --restart-step 0.")
    if p_e is not None and rho_g_of_p is not None:
        pg_e = p_e + (md.restart["PCOW"][gi] if (args.capillary == "deck" and "PCOW" in (md.restart or {})) else 0.0)
        giip_sm3 = float((m.elem_vars["pore_volume"] * sg_e * rho_g_of_p(pg_e)).sum() / rho_g_sc)
        unit = eclio.UNITS[md.units]["gas_rate"]
        rep("- Initial gas in place in the MOOSE model: %.6g sm3 = %.6g %s (surface conditions, from PVDG)."
            % (giip_sm3, giip_sm3 / unit, "Mscf" if md.units == "FIELD" else "sm3"))
    if p_e is not None:
        m.elem_vars["p_init"] = p_e
        m.elem_vars["sg_init"] = sg_e
        m.node_vars["p_init_nodal"] = p_n
        m.node_vars["sg_init_nodal"] = np.clip(sg_n, 0.0, 0.999)

    # --- Pc = 0 surrogate: water held above residual by capillarity in Eclipse is made immobile in MOOSE
    #     (gas-cap cells binned by initial Sw into sub-blocks whose water curve starts at that Sw;
    #      gas curves unchanged) -> static initial state, gas in place preserved exactly.
    if sg_e is not None and not args.mobile_connate_water and init_mode == "restart" and args.capillary == "none":
        base = {r_["block"]: r_ for r_ in relperm}
        names_by_id = dict(m.block_names)
        new_blocks = m.blocks.copy()
        new_names = dict(names_by_id)
        extra = []
        nid = max(names_by_id) + 1 if names_by_id else 1
        sw_e = 1.0 - sg_e
        combos = {}
        for b_ in sorted(set(m.blocks.tolist())):
            rp_ = base.get(names_by_id.get(b_))
            if rp_ is None:
                continue
            sel = (m.blocks == b_) & (sg_e > 1e-3) & (sw_e > rp_["swr"] + 1e-3)
            if not sel.any():
                continue
            width = args.swi_bin
            upper = rp_["swr"] + np.ceil((sw_e[sel] - rp_["swr"]) / width) * width
            upper = np.minimum(upper, 0.99 - rp_["sgr"])
            idx = np.nonzero(sel)[0]
            upper = np.round(upper, 3)
            for u in sorted(set(upper.tolist())):
                ids = idx[upper == u]
                nm = "%s_swi%04d" % (names_by_id[b_], int(round(u * 1000)))
                new_blocks[ids] = nid
                new_names[nid] = nm
                r2 = dict(rp_)
                r2.update(block=nm, swr_water=u)
                extra.append(r2)
                nid += 1
        if extra:
            m.blocks = new_blocks
            used = set(m.blocks.tolist())
            m.block_names = {b_: n_ for b_, n_ in new_names.items() if b_ in used}
            relperm = [r_ for r_ in relperm if r_["block"] in m.block_names.values()] + extra
            n_cells = int(sum((m.blocks == b_).sum() for b_, n_ in m.block_names.items() if "_swi" in n_))
            rep("- Pc = 0 surrogate: %d gas-zone cells hold water above the table residual (Eclipse keeps it there with "
                "capillary pressure). They are put in %d sub-blocks (initial Sw binned by %.3f) whose water curve "
                "starts at that Sw, so this water is immobile until it is exceeded; gas curves are unchanged. "
                "Use --mobile-connate-water to switch this off." % (n_cells, len(extra), args.swi_bin))

    # --- capillary pressure from the deck (Brooks-Corey fit per table, Eclipse end-point scaling by
    #     sub-blocks binned on the cell's PCW (Pc at SWLPC) and SWLPC)
    if args.capillary == "deck":
        relperm, cap_lines = _capillary_blocks(m, relperm, deck, P, gi, sg_e, md.restart, args)
        for line in cap_lines:
            rep(line)

    # --- aquifers -> sidesets + (finite) Fetkovich model
    aquifers = []
    if deck:
        Pu = eclio.UNITS[deck.units]
        for aid, aq in sorted(deck.aquifers().items()):
            for cnx in aq["connections"]:
                face = cnx["face"].replace("X", "I").replace("Y", "J").replace("Z", "K")
                msk = np.zeros((grid.nz, grid.ny, grid.nx), bool)
                k1_, k2_ = cnx["k1"], cnx["k2"]
                if md.kmap is not None:
                    k1_, k2_ = md.kmap[min(k1_, len(md.kmap) - 1)], md.kmap[min(k2_, len(md.kmap) - 1)]
                msk[k1_ - 1:k2_, cnx["j1"] - 1:cnx["j2"], cnx["i1"] - 1:cnx["i2"]] = True
                nm = "aquifer_%d_%s" % (aid, face.replace("-", "m").replace("+", "p"))
                prev = m.sidesets.get(nm)
                add_face_sideset(grid, m, nm, msk.ravel(), face)
                if prev is not None:  # merge repeated faces for same aquifer
                    e = np.concatenate([prev[0], m.sidesets[nm][0]])
                    s_ = np.concatenate([prev[1], m.sidesets[nm][1]])
                    u, ix = np.unique(np.stack([e, s_], 1), axis=0, return_index=True)
                    m.sidesets[nm] = (e[np.sort(ix)], s_[np.sort(ix)])
            names = sorted({k for k in m.sidesets if k.startswith("aquifer_%d_" % aid)})
            names = [k for k in names if len(m.sidesets[k][0])]
            area = sum(sideset_area(m, k) for k in names)
            if not names or area <= 0:
                rep("- Aquifer %d: no external faces found for its AQUANCON connections; skipped." % aid)
                continue
            el = np.unique(np.concatenate([m.sidesets[k][0] for k in names]))
            aq.update(id=aid, boundaries=names, area=area)
            if not aq.get("datum_depth") or aq["datum_depth"] <= 0:
                aq["datum_depth"] = float(depth_e[el].mean())
                rep("- Aquifer %d: datum depth defaulted -> mean depth of connected cells %.1f m." % (aid, aq["datum_depth"]))
            if aq.get("p0") is None:
                if p_e is not None:
                    elw = el[sg_e[el] < 0.01] if (sg_e[el] < 0.01).any() else el   # water-zone cells only
                    aq["p0"] = float(np.mean(p_e[elw] + rho_w * G * (aq["datum_depth"] - depth_e[elw])))
                    rep("- Aquifer %d: initial pressure defaulted -> equilibrium with connected cells %.4g Pa." % (aid, aq["p0"]))
                else:
                    aq["p0"] = 1.5e7
                    rep("- Aquifer %d: initial pressure unknown -> 15 MPa placeholder. EDIT." % aid)
            mult = None
            if aq["type"] == "FETKOVICH" and isinstance(aq.get("pi_raw"), float):
                pi_si = aq["pi_raw"] * Pu["liq_rate"] / 86400.0 / Pu["pressure"]  # m3/s/Pa
                mult = pi_si * rho_w_sc / area
                how = "Fetkovich PI %.4g (deck units) spread over %.4g m2" % (aq["pi_raw"], area)
                if aq.get("volume") and aq.get("ct"):
                    aq["capacity"] = aq["volume"] * aq["ct"] * rho_w_sc   # kg/Pa
                    how += "; finite: p_aq drops by cumulative influx / (rho V ct), V ct = %.4g m3/Pa" % (aq["volume"] * aq["ct"])
            elif aq["type"] == "CARTER_TRACY" and isinstance(aq.get("perm_md"), float) and aq.get("r0"):
                mult = rho_w * aq["perm_md"] * eclio.MD_TO_M2 / (args.mu_w * aq["r0"])
                how = "Carter-Tracy approximated as steady flux rho*k/(mu*r0), constant pressure"
            if mult is None:
                mult = args.aquifer_mult
                how = "no usable aquifer data: default multiplier %g" % mult
            aq["mult"] = mult
            aquifers.append(aq)
            rep("- Aquifer %d (%s): boundaries %s; flux = %.3g kg/s/m2/Pa x (p - p_aq(z,t)); %s."
                % (aid, aq["type"], " ".join(names), mult, how))

    # --- wells
    wells = []
    if deck:
        sch = deck.schedule()
        lut = -np.ones(n, np.int64)
        lut[gi] = np.arange(len(gi))
        cen = cell_centres(m)
        rows = []
        for wname, cons in sch["compdat"].items():
            conns = []
            for cc in cons:
                if cc["status"] not in ("OPEN", "AUTO"):
                    continue
                kk_ = md.kmap[cc["k"]] if md.kmap is not None else cc["k"]
                g = grid.gidx(cc["i"] - 1, cc["j"] - 1, kk_ - 1)
                e = lut[g]
                if e < 0:
                    rep("- Well %s connection (%d,%d,%d) is in an inactive/dropped cell; ignored."
                        % (wname, cc["i"], cc["j"], cc["k"]))
                    continue
                if isinstance(cc.get("kh"), float) and cc["kh"] > 0:
                    wgt = cc["kh"]
                elif isinstance(cc.get("cf"), float) and cc["cf"] > 0:
                    wgt = cc["cf"]
                else:
                    h = (m.nodes[m.conn[e, 4:]].mean(0) - m.nodes[m.conn[e, :4]].mean(0))[2]
                    wgt = abs(h) * math.sqrt(kx[g] * ky[g])
                cf_ = cc["cf"] if isinstance(cc.get("cf"), float) and cc["cf"] > 0 else None
                dup = [c_ for c_ in conns if c_["elem"] == int(e)]
                if dup:
                    dup[0]["w"] += wgt
                    dup[0]["cf"] = None if (cf_ is None or dup[0]["cf"] is None) else dup[0]["cf"] + cf_
                else:
                    conns.append(dict(cell=(cc["i"], cc["j"], kk_), elem=int(e), xyz=cen[e], w=wgt, cf=cf_))
            ev = sch["events"].get(wname, [])
            if not conns:
                if ev:
                    rep("- Well %s has rate data but no usable connections; skipped." % wname)
                continue
            tw = sum(c_["w"] for c_ in conns) or 1.0
            for c_ in conns:
                c_["frac"] = c_["w"] / tw
                rows.append([wname, *c_["cell"], *(c_["xyz"] + m.origin), c_["frac"]])
            # rate history -> mass rates kg/s (+ in, - out)
            ts, qg, qw = [], [], []
            for (t, kind, g_, w_, st, src) in ev:
                sgn = 1.0 if kind == "INJ" else -1.0
                if ts and abs(ts[-1] - t) < 1e-6:
                    ts.pop(); qg.pop(); qw.pop()
                ts.append(t)
                qg.append(sgn * g_ * rho_g_sc / 86400.0)
                qw.append(sgn * w_ * rho_w_sc / 86400.0)
            if not ts or (not any(qg) and not any(qw)):
                rep("- Well %s: no non-zero rates in schedule; connections written but no source term." % wname)
            if ts and ts[0] > 0:   # no flow before the first schedule entry (PiecewiseConstant holds y[0] before x[0])
                ts.insert(0, 0.0); qg.insert(0, 0.0); qw.insert(0, 0.0)
            if args.average_rates and len(ts) > 1:
                # mass-conserving averaging of the rate history over periods of N days
                T_end = sch["end"]
                per = args.average_rates * 86400.0
                edges = np.arange(0.0, T_end + per, per)
                edges[-1] = max(edges[-1], T_end)
                tt = np.array(ts + [max(T_end, ts[-1])])

                def cum_at(t_, q_):
                    q_ = np.array(q_)
                    seg = np.clip(np.minimum(tt[1:], t_) - tt[:-1], 0, None)
                    return float((q_ * seg).sum())
                new_t, new_g, new_w = [], [], []
                for a_, b_ in zip(edges[:-1], edges[1:]):
                    new_t.append(a_)
                    new_g.append((cum_at(b_, qg) - cum_at(a_, qg)) / (b_ - a_))
                    new_w.append((cum_at(b_, qw) - cum_at(a_, qw)) / (b_ - a_))
                # drop repeated values
                keep = [0] + [k_ for k_ in range(1, len(new_t)) if (new_g[k_], new_w[k_]) != (new_g[k_ - 1], new_w[k_ - 1])]
                ts = [new_t[k_] for k_ in keep]; qg = [new_g[k_] for k_ in keep]; qw = [new_w[k_] for k_ in keep]
            if len(ts) == 1:   # PiecewiseConstant needs two points
                ts.append(max(sch["end"], ts[0] + 86400.0)); qg.append(qg[0]); qw.append(qw[0])
            ref_ft = (sch.get("welspecs", {}).get(wname) or {}).get("ref_depth")
            Lu_ = 0.3048 if deck.units == "FIELD" else 1.0
            wells.append(dict(name=wname, conns=conns, t=ts, qg=qg, qw=qw,
                              z_ref=(-float(ref_ft) * Lu_) if ref_ft else float(np.mean([c_["xyz"][2] for c_ in conns]))))
        end_time = sch["end"]
        if wells:
            rep("- Wells: %s. Rates from WCONHIST/WCONPROD/WCONINJE/WCONINJH converted to mass rates with "
                "surface densities gas %.3f, water %.1f kg/m3; split over connections by Kh/CF."
                % (", ".join(w["name"] for w in wells), rho_g_sc, rho_w_sc))
            if args.average_rates:
                rep("- Well rates averaged over %g-day periods (cumulative volumes preserved at period ends)." % args.average_rates)
            if args.fv and args.well_model == "connection":
                rep("- Wells (gas production): Eclipse-style connections. Each completion cell gets the share "
                    "CF_c rho krg/mu_g / sum(CF rho krg/mu_g) of the target rate (mobility from the last time step), "
                    "limited to what it can deliver with the bottom-hole pressure at %.4g Pa (1 atm, WCONHIST default) "
                    "at the BHP reference depth: q_c <= CF_c rho krg/mu_g (p_c - p_bhp(z_c)). A completion whose gas is "
                    "trapped or gone produces nothing; if no completion can deliver, the well produces less than the "
                    "history (reported as <well>_gas_rate_actual_kgps). WBP is the CF-weighted gas pressure of the "
                    "completion cells corrected to the BHP reference depth with the gas density (WPAVE ... RES)." % 101325.0)
                if getattr(args, "well_water", "on") == "on":
                    rep("- Water production (--well-water on): each completion also produces water at the same drawdown "
                        "as its gas, CF_c rho_w krw/mu_w x dp_c (water mobility from the last time step), as Eclipse "
                        "produces every mobile phase. The history gas rate stays the target (Eclipse controlled RESV, "
                        "reservoir volume of gas + water; the difference is the water volume, < 0.5 %% for Croft). "
                        "Reported as <well>_water_rate_actual_kgps and <well>_cum_water_actual_kg; compare with Eclipse WWPR.")
                rep("  NOTE: injection and history water rates are still rate-controlled point sources.")
            else:
                rep("  NOTE: wells are rate-controlled point sources for the controlled phase (gas for GRAT); BHP limits are "
                    "NOT enforced and water is only produced where a water rate is given.")
        with open(args.output + "_wells.csv", "w") as f:
            f.write("well,i,j,k,x,y,z_elev,rate_fraction\n")
            for r in rows:
                f.write("%s,%d,%d,%d,%.3f,%.3f,%.3f,%.6f\n" % tuple(r))
    else:
        end_time = args.end_time * 86400.0

    # --- per-well source weights for the finite-volume input (rate fraction / cell volume, 1/m3)
    for wl in wells:
        arr = np.zeros(len(m.conn))
        for c_ in wl["conns"]:
            arr[c_["elem"]] += c_["frac"] / vol[c_["elem"]]
        wl["weight_var"] = ("ww_" + safe(wl["name"]))[:32]
        m.elem_vars[wl["weight_var"]] = arr
        if all(c_.get("cf") for c_ in wl["conns"]):
            cfu = (1e-3 * 0.158987294928 / 86400.0 / 6894.757293168) if deck.units == "FIELD" else (1e-3 / 86400.0 / 1e5)
            arr2 = np.zeros(len(m.conn))
            for c_ in wl["conns"]:
                arr2[c_["elem"]] += c_["cf"] * cfu / vol[c_["elem"]]
            wl["cf_var"] = ("wc_" + safe(wl["name"]))[:32]
            m.elem_vars[wl["cf_var"]] = arr2
            wl["cf_total"] = float(sum(c_["cf"] for c_ in wl["conns"]) * cfu)
    if "sg_init" in m.elem_vars:
        m.elem_vars["sg_init"] = np.clip(m.elem_vars["sg_init"], 0.0, 0.999)

    # --- write mesh
    log("Writing %s.e ..." % args.output)
    perm_exo = write_exodus(args.output + ".e", m, title="petrel2moose: %s" % os.path.basename(args.case or args.egrid or args.grdecl))
    nnc_file = None
    if use_nnc:
        A, B, Tn, st = eclipse_nncs(args.egrid, args.init, md.units, md.kmap, m, log)
        if len(A):
            inv_exo = np.empty(len(perm_exo), np.int64)
            inv_exo[perm_exo] = np.arange(len(perm_exo))
            nnc_file = args.output + "_nnc.csv"
            with open(nnc_file, "w") as f:
                f.write("elem_a,elem_b,trans,i_a,j_a,k_a,i_b,j_b,k_b\n")
                for a_, b_, t_ in zip(A, B, Tn):
                    ca, cb = m.cells[a_] + 1, m.cells[b_] + 1
                    f.write("%d,%d,%.9e,%d,%d,%d,%d,%d,%d\n" % (inv_exo[a_], inv_exo[b_], t_ * args.nnc_mult,
                                                              ca[0], ca[1], ca[2], cb[0], cb[1], cb[2]))
            rep("- Non-neighbour connections: %d MOOSE element pairs that share no mesh face (fault throw, "
                "pinch-outs) are connected with Eclipse's own transmissibilities (TRANX/Y/Z and NNC/TRANNNC: %d "
                "Eclipse connections, %d of them NNCs; %.2f%% of the total transmissibility%s), written to %s "
                "and applied by FVPorousFlowNNCFlux (element ids in the Exodus file order; Mesh "
                "allow_renumbering = false)."
                % (len(A), st["n_fine"], st["n_eclipse_nnc"], 100 * st["frac_trans"],
                   "" if md.kmap is None else ", summed over merged layers", os.path.basename(nnc_file),))

    # --- write MOOSE input
    rock = None
    if deck and deck.kw.get("ROCK") and not args.no_rock_compressibility:
        r0 = deck.kw["ROCK"][0][0]
        Pu_ = eclio.UNITS[deck.units]["pressure"]
        if len(r0) >= 2 and r0[1] and r0[1] > 0:
            rock = dict(p_ref=float(r0[0]) * Pu_, c_r=float(r0[1]) / Pu_)
            rep("- Rock compressibility (ROCK): %.4g 1/Pa about %.4g Pa%s." % (
                rock["c_r"], rock["p_ref"], "" if args.fv else " (used by the --fv input only)"))
    params = dict(T=T, rho_w=rho_w, pvtw=pvtw, rho_w_sc=rho_w_sc, rock=rock,
                  gas_pvt_file=os.path.basename(gas_pvt_file) if gas_pvt_file else None, rho_g_res=rho_g_res, swr=swr, relperm=relperm, capillary=args.capillary, eq=eq, init_mode=init_mode, aquifers=aquifers, wells=wells,
                  end_time=end_time, mesh=os.path.basename(args.output) + ".e", origin=m.origin,
                  zmin=float(m.nodes[:, 2].min()), zmax=float(m.nodes[:, 2].max()),
                  sch_times=(deck.schedule()["times"] if deck else []),
                  nnc_file=os.path.basename(nnc_file) if nnc_file else None)
    with open(args.output + ".i", "w") as f:
        f.write(moose_input_fv(params, args) if args.fv else moose_input(params, args))

    with open(args.output + "_report.md", "w") as f:
        f.write("# petrel2moose conversion report\n\n")
        f.write("Source: `%s`\n\n" % (args.case or args.egrid or args.grdecl))
        f.write("\n".join(report) + "\n")
    with open(args.output + "_meta.json", "w") as f:
        json.dump(dict(origin=m.origin.tolist(), stats=m.stats, units_in=md.units,
                       sidesets={k: int(len(v[0])) for k, v in m.sidesets.items()},
                       blocks={str(k): v for k, v in m.block_names.items()},
                       rho_water_sc=float(params["rho_w_sc"]), rho_gas_sc=float(locals().get("rho_g_sc") or params.get("rho_g_sc", 0) or 0) or None), f, indent=2)
    log("\n".join(report))
    return m, params



HEX_FACES = np.array([[0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7], [0, 3, 2, 1], [4, 5, 6, 7]])


def mesh_face_pairs(conn):
    """Set of element pairs (a < b) that share a whole face in the mesh."""
    f = np.sort(conn[:, HEX_FACES].reshape(-1, 4), axis=1)
    owner = np.repeat(np.arange(len(conn)), 6)
    _, inv, cnt = np.unique(f, axis=0, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    sel = cnt[inv] == 2
    order = np.argsort(inv[sel], kind="stable")
    o = owner[sel][order].reshape(-1, 2)
    a, b = np.minimum(o[:, 0], o[:, 1]), np.maximum(o[:, 0], o[:, 1])
    return set(zip(a.tolist(), b.tolist()))


def eclipse_nncs(egrid, init, units, kmap, m, log=print):
    """Eclipse connections (TRANX/Y/Z and NNC/TRANNNC, as Eclipse computed them) between pairs of MOOSE
    elements that do not share a mesh face, e.g. across faults with throw, where the MOOSE mesh has no
    face and would otherwise give no flow. Returns (elem_a, elem_b, trans_m3) with element indices in
    mesh order (merged layers: transmissibilities of all fine connections between the same pair summed)."""
    items = eclio.read_unformatted(egrid)
    gh = eclio.first(items, "GRIDHEAD")
    nx, ny, nz = int(gh[1]), int(gh[2]), int(gh[3])
    act = np.asarray(eclio.first(items, "ACTNUM")).ravel() > 0
    gact = np.nonzero(act)[0]
    it = eclio.read_unformatted(init)
    u = eclio.UNITS[units]
    tsi = 1.0e-3 * u["liq_rate"] / 86400.0 / u["pressure"]   # Eclipse T -> m^3 (Pa s * m^3/s / Pa)
    conns = []
    for kw, step in (("TRANX", 1), ("TRANY", nx), ("TRANZ", nx * ny)):
        t = eclio.first(it, kw)
        if t is None:
            continue
        t = np.asarray(t, float)
        g1 = gact
        g2 = g1 + step
        i, j, k = g1 % nx, (g1 // nx) % ny, g1 // (nx * ny)
        ok = t > 0
        ok &= {1: i < nx - 1, nx: j < ny - 1, nx * ny: k < nz - 1}[step]
        g2c = np.clip(g2, 0, len(act) - 1)
        ok &= act[g2c]
        conns.append((g1[ok], g2[ok], t[ok]))
    n1, n2, tn = eclio.first(items, "NNC1"), eclio.first(items, "NNC2"), eclio.first(it, "TRANNNC")
    n_nnc = 0
    if n1 is not None and tn is not None and len(n1) == len(tn):
        n1, n2, tn = np.asarray(n1, int) - 1, np.asarray(n2, int) - 1, np.asarray(tn, float)
        ok = tn > 0
        conns.append((n1[ok], n2[ok], tn[ok]))
        n_nnc = int(ok.sum())
    g1 = np.concatenate([c[0] for c in conns])
    g2 = np.concatenate([c[1] for c in conns])
    T = np.concatenate([c[2] for c in conns]) * tsi

    # fine global cell -> MOOSE element (layers merged through kmap)
    nzc = int(kmap.max()) if kmap is not None else nz
    lookup = np.full(nx * ny * nzc, -1, np.int64)
    lookup[np.asarray(m.gidx, np.int64)] = np.arange(len(m.gidx))

    def to_elem(g):
        i, j, k = g % nx, (g // nx) % ny, g // (nx * ny)
        kc = (kmap[k + 1] - 1) if kmap is not None else k
        return lookup[kc * nx * ny + j * nx + i]

    e1, e2 = to_elem(g1), to_elem(g2)
    ok = (e1 >= 0) & (e2 >= 0) & (e1 != e2)
    e1, e2, T = e1[ok], e2[ok], T[ok]
    a, b = np.minimum(e1, e2), np.maximum(e1, e2)
    nb = mesh_face_pairs(m.conn)
    key = a * len(m.gidx) + b
    isnb = np.fromiter(((x, y) in nb for x, y in zip(a.tolist(), b.tolist())), bool, len(a))
    t_all = T.sum()
    a, b, T, key = a[~isnb], b[~isnb], T[~isnb], key[~isnb]
    uk, inv = np.unique(key, return_inverse=True)
    Ts = np.bincount(inv.ravel(), weights=T)
    A, B = uk // len(m.gidx), uk % len(m.gidx)
    stats = dict(n_fine=int(len(g1)), n_eclipse_nnc=n_nnc, n_pairs=int(len(uk)),
                 frac_trans=float(Ts.sum() / max(t_all, 1e-300)))
    log("  NNC: %d Eclipse connections (%d from NNC/TRANNNC); %d element pairs without a mesh face get a "
        "non-neighbour connection (%.2f%% of total transmissibility)"
        % (stats["n_fine"], n_nnc, stats["n_pairs"], 100 * stats["frac_trans"]))
    return A, B, Ts, stats


def _eclipse_pv(md, phi, keep):
    """Active pore volume of the original grid (geometric), for the 'kept' check."""
    act = (md.grid.actnum.ravel() > 0) & (phi > 0)
    g = np.nonzero(act)[0]
    k, rem = np.divmod(g, md.grid.nx * md.grid.ny)
    j, i = np.divmod(rem, md.grid.nx)
    Pc, _ = md.grid.cell_corners(np.stack([i, j, k], 1))
    return float((np.abs(_hex_volume(Pc)) * phi[g]).sum())


def _hex_volume(P):
    from cornerpoint import GAUSS8, hex_jacobians
    return hex_jacobians(P, GAUSS8).sum(1)


# ----------------------------------------------------------------------------------------------
# MOOSE input file
# ----------------------------------------------------------------------------------------------
def moose_input(p, args):
    L = []
    w = L.append
    eq = p["eq"]
    w("# =============================================================================================")
    w("# MOOSE PorousFlow model generated by petrel2moose.py")
    w("#   2-phase immiscible: phase 0 = water (Water97), phase 1 = gas (methane), isothermal")
    w("#   primary variables: pwater (water porepressure, Pa) and sgas (gas saturation)")
    w("#   mesh: %s   (x,y shifted by origin %.1f %.1f;  z = elevation, m)" % (p["mesh"], p["origin"][0], p["origin"][1]))
    w("#   Run:  porous_flow-opt -i <this file>     (or any MOOSE app that includes PorousFlow)")
    w("# =============================================================================================")
    w("")
    w("[Mesh]")
    w("  [file]")
    w("    type = FileMeshGenerator")
    w("    file = %s" % p["mesh"])
    w("    use_for_exodus_restart = true")
    w("  []")
    w("[]")
    w("")
    w("[GlobalParams]")
    w("  PorousFlowDictator = dictator")
    w("  gravity = '0 0 -9.81'")
    w("[]")
    w("")
    w("[UserObjects]")
    w("  [dictator]")
    w("    type = PorousFlowDictator")
    w("    porous_flow_vars = 'pwater sgas'")
    w("    number_fluid_phases = 2")
    w("    number_fluid_components = 2")
    w("  []")
    cap = p.get("capillary") == "deck" and any("pc" in rp for rp in p["relperm"])
    if not cap:
        w("  [pc]")
        w("    type = PorousFlowCapillaryPressureConst")
        w("    pc = 0")
        w("  []")
    else:
        for rp in p["relperm"]:
            w("  [pc_%s]" % rp["block"])
            if "pc" in rp:
                c = rp["pc"]
                w("    type = PorousFlowCapillaryPressureBC # Brooks-Corey fit to SATNUM table %d, Pc(Swl) %s Pa"
                  % (rp["table"], fmt(round(c["pcw"], 1))))
                w("    pe = %s" % fmt(round(c["pe"], 2)))
                w("    lambda = %s" % fmt(round(c["lam"], 4)))
                w("    sat_lr = %s" % fmt(round(c["sat_lr"], 4)))
                w("    pc_max = %s" % fmt(round(c["pc_max"], 1)))
            else:
                w("    type = PorousFlowCapillaryPressureConst")
                w("    pc = 0")
            w("  []")
    w("[]")
    w("")
    # Variables + ICs
    w("[Variables]")
    for v, f in (("pwater", "p_init_nodal"), ("sgas", "sg_init_nodal")):
        w("  [%s]" % v)
        if p["init_mode"] == "restart":
            w("    initial_from_file_var = %s" % f)
        w("  []")
    w("[]")
    w("")
    fn, ic = [], []   # [Functions] and [ICs] collected, written once each
    if p["init_mode"] != "restart":
        if eq:
            gwc_z, pg = -eq["gwc"], eq["p_gwc"]
            fn.append("  # hydrostatic initialisation from EQUIL: gas column above GWC, water below, zero Pc")
        else:
            gwc_z, pg = 0.5 * (p["zmin"] + p["zmax"]), 1.5e7
            fn.append("  # !!! No EQUIL / restart available: EDIT gwc_elev and p_gwc below before running !!!")
        fn += ["  [p_ic]", "    type = ParsedFunction", "    symbol_names = 'gwc_elev p_gwc rho_w rho_g g'",
               "    symbol_values = '%s %s %s %s 9.81'" % (fmt(gwc_z), fmt(pg), fmt(p["rho_w"]), fmt(p["rho_g_res"] or 150.0)),
               "    expression = 'if(z > gwc_elev, p_gwc - rho_g * g * (z - gwc_elev), p_gwc - rho_w * g * (z - gwc_elev))'",
               "  []",
               "  [sg_ic]", "    type = ParsedFunction", "    symbol_names = 'gwc_elev sg_max'",
               "    symbol_values = '%s %s'" % (fmt(gwc_z), fmt(1.0 - p["swr"] - 1e-3)),
               "    expression = 'if(z > gwc_elev, sg_max, 0)'", "  []"]
        ic += ["  [p]", "    type = FunctionIC", "    variable = pwater", "    function = p_ic", "  []",
               "  [sg]", "    type = FunctionIC", "    variable = sgas", "    function = sg_ic", "  []"]
    for aq in p["aquifers"]:
        a = aq["id"]
        if aq.get("capacity"):
            fn += ["  [p_aq%d_fn]   # finite (Fetkovich) aquifer %d: hydrostatic about datum, plus cumulative "
                   "outflow / (rho V ct)" % (a, a),
                   "    type = ParsedFunction", "    symbol_names = 'p0 datum_elev rho_w g cum cap'",
                   "    symbol_values = '%s %s %s 9.81 aq%d_cum_out_kg %s'"
                   % (fmt(aq["p0"]), fmt(-aq["datum_depth"]), fmt(p["rho_w"]), a, fmt(aq["capacity"])),
                   "    expression = 'p0 - rho_w * g * (z - datum_elev) + cum / cap'", "  []"]
        else:
            fn += ["  [p_aq%d_fn]   # aquifer %d pressure, hydrostatic about its datum (constant in time)" % (a, a),
                   "    type = ParsedFunction", "    symbol_names = 'p0 datum_elev rho_w g'",
                   "    symbol_values = '%s %s %s 9.81'" % (fmt(aq["p0"]), fmt(-aq["datum_depth"]), fmt(p["rho_w"])),
                   "    expression = 'p0 - rho_w * g * (z - datum_elev)'", "  []"]
    w_functions_wells(fn.append, p)
    if fn:
        w("[Functions]")
        L.extend(fn)
        w("[]")
        w("")
    if ic:
        w("[ICs]")
        L.extend(ic)
        w("[]")
        w("")

    w("[AuxVariables]")
    for nm in ("poro", "perm_xx", "perm_yy", "perm_zz", "perm_xy"):
        w("  [%s]" % nm)
        w("    family = MONOMIAL")
        w("    order = CONSTANT")
        w("    initial_from_file_var = %s" % nm)
        w("  []")
    w("  [mass_frac_phase0_species0]")
    w("    initial_condition = 1 # water phase is pure water")
    w("  []")
    w("  [mass_frac_phase1_species0]")
    w("    initial_condition = 0 # gas phase is pure methane")
    w("  []")
    for aq in p["aquifers"]:
        w("  [p_aq%d]   # aquifer %d pressure seen by its boundary [Pa]" % (aq["id"], aq["id"]))
        w("  []")
        w("  [aq%d_flux]   # nodal mass outflow to aquifer %d [kg/s] (save_in)" % (aq["id"], aq["id"]))
        w("  []")
    for nm in ("swater", "gas_density"):
        w("  [%s]" % nm)
        w("    family = MONOMIAL")
        w("    order = CONSTANT")
        w("  []")
    w("[]")
    w("")
    w("[AuxKernels]")
    for aq in p["aquifers"]:
        w("  [p_aq%d]" % aq["id"])
        w("    type = FunctionAux")
        w("    variable = p_aq%d" % aq["id"])
        w("    function = p_aq%d_fn" % aq["id"])
        w("    execute_on = 'initial timestep_begin'")
        w("  []")
    for nm, prop, ph in (("swater", "saturation", 0), ("gas_density", "density", 1)):
        w("  [%s]" % nm)
        w("    type = PorousFlowPropertyAux")
        w("    variable = %s" % nm)
        w("    property = %s" % prop)
        w("    phase = %d" % ph)
        w("    execute_on = timestep_end")
        w("  []")
    w("[]")
    w("")

    w("[Kernels]")
    w("  [mass_water_dot]")
    w("    type = PorousFlowMassTimeDerivative")
    w("    fluid_component = 0")
    w("    variable = pwater")
    w("  []")
    w("  [flux_water]")
    w("    type = PorousFlowAdvectiveFlux")
    if cap:
        w("    full_upwind_threshold = 100 # keep full upwinding (no switch to the non-conservative fallback)")
    w("    fluid_component = 0")
    w("    variable = pwater")
    w("  []")
    w("  [mass_gas_dot]")
    w("    type = PorousFlowMassTimeDerivative")
    w("    fluid_component = 1")
    w("    variable = sgas")
    w("  []")
    w("  [flux_gas]")
    w("    type = PorousFlowAdvectiveFlux")
    if cap:
        w("    full_upwind_threshold = 100 # keep full upwinding (no switch to the non-conservative fallback)")
    w("    fluid_component = 1")
    w("    variable = sgas")
    w("  []")
    w("[]")
    w("")
    w("[FluidProperties]")
    if p.get("pvtw"):
        pw = p["pvtw"]
        rho_ref = p["rho_w_sc"] / pw["bw"]
        w("  [water]   # from PVTW: rho = rho_sc/Bw * exp(cw (p - pref)), constant viscosity")
        w("    type = SimpleFluidProperties")
        w("    density0 = %s" % fmt(rho_ref * math.exp(-pw["cw"] * pw["pref"])))
        w("    bulk_modulus = %s" % fmt(1.0 / pw["cw"]))
        w("    viscosity = %s" % fmt(pw["mu"]))
        w("    thermal_expansion = 0")
        w("  []")
    else:
        w("  [water]")
        w("    type = Water97FluidProperties")
        w("  []")
    if p.get("gas_pvt_file"):
        w("  [gas]   # from PVDG (Eclipse dry-gas table)")
        w("    type = TabulatedBicubicFluidProperties")
        w("    fluid_property_file = %s" % p["gas_pvt_file"])
        w("    interpolated_properties = 'density viscosity'")
        w("    out_of_bounds_behavior = set_to_closest_bound")
        w("  []")
    else:
        w("  [gas]")
        w("    type = MethaneFluidProperties")
        w("  []")
    w("[]")
    w("")
    w("[Materials]")
    w("  [temperature]")
    w("    type = PorousFlowTemperature")
    w("    temperature = %s # K" % fmt(p["T"]))
    w("  []")
    if not cap:
        w("  [ppss]")
        w("    type = PorousFlow2PhasePS")
        w("    phase0_porepressure = pwater")
        w("    phase1_saturation = sgas")
        w("    capillary_pressure = pc")
        w("  []")
    else:
        for rp in p["relperm"]:
            w("  [ppss_%s]" % rp["block"])
            w("    type = PorousFlow2PhasePS")
            w("    block = %s" % rp["block"])
            w("    phase0_porepressure = pwater")
            w("    phase1_saturation = sgas")
            w("    capillary_pressure = pc_%s" % rp["block"])
            w("  []")
    w("  [massfrac]")
    w("    type = PorousFlowMassFraction")
    w("    mass_fraction_vars = 'mass_frac_phase0_species0 mass_frac_phase1_species0'")
    w("  []")
    for ph, nm in ((0, "water"), (1, "gas")):
        w("  [%s]" % nm)
        w("    type = PorousFlowSingleComponentFluid")
        w("    fp = %s" % nm)
        w("    phase = %d" % ph)
        w("    compute_enthalpy = false")
        w("    compute_internal_energy = false")
        w("  []")
    for rp in p["relperm"]:
        blk = rp["block"]
        multi = len(p["relperm"]) > 1
        for ph, nm in ((0, "water"), (1, "gas")):
            w("  [relperm_%s%s]" % (nm, "_" + blk if multi else ""))
            w("    type = PorousFlowRelativePermeabilityCorey")
            if multi:
                w("    block = %s" % blk)
            w("    phase = %d" % ph)
            w("    n = %s" % fmt(round(rp["nw"] if ph == 0 else rp["ng"], 3)))
            sw_res = rp.get("swr_water", rp["swr"])
            w("    s_res = %s" % fmt(round(sw_res if ph == 0 else rp["sgr"], 4)))
            w("    sum_s_res = %s" % fmt(round((sw_res if ph == 0 else rp.get("swr_gas", rp["swr"])) + rp["sgr"], 4)))
            w("    scaling = %s # end point (SATNUM table %d)" % (fmt(round(rp["kw_max"] if ph == 0 else rp["kg_max"], 4)), rp["table"]))
            w("  []")
    w("  [porosity]")
    w("    type = PorousFlowPorosityConst")
    w("    porosity = poro")
    w("  []")
    w("  [permeability]")
    w("    type = PorousFlowPermeabilityConstFromVar")
    w("    perm_xx = perm_xx")
    w("    perm_yy = perm_yy")
    w("    perm_zz = perm_zz")
    w("    perm_xy = perm_xy")
    w("    perm_yx = perm_xy")
    w("  []")
    w("[]")
    w("")
    if p["aquifers"]:
        w("[BCs]")
        for aq in p["aquifers"]:
            mlt = aq["mult"]
            w("  # aquifer %d (%s): flux out = %s * (pwater - p_aq%d)  [kg/s/m2]; negative = influx"
              % (aq["id"], aq["type"], fmt(mlt), aq["id"]))
            w("  [aquifer_%d]" % aq["id"])
            w("    type = PorousFlowPiecewiseLinearSink")
            w("    variable = pwater")
            w("    boundary = '%s'" % " ".join(aq["boundaries"]))
            w("    fluid_phase = 0")
            w("    PT_shift = p_aq%d" % aq["id"])
            w("    pt_vals = '-1e9 1e9'")
            w("    multipliers = '%s %s'" % (fmt(-1e9 * mlt), fmt(1e9 * mlt)))
            w("    mass_fraction_component = 0")
            w("    save_in = aq%d_flux" % aq["id"])
            w("  []")
        w("[]")
        w("")
    if p["wells"]:
        w("[DiracKernels]")
        for wl in p["wells"]:
            nm = safe(wl["name"])
            for ci, c in enumerate(wl["conns"]):
                for ph, var, q in (("gas", "sgas", wl["qg"]), ("water", "pwater", wl["qw"])):
                    if not any(q):
                        continue
                    w("  [%s_%s_c%d]   # cell %d,%d,%d" % (nm, ph, ci, *c["cell"]))
                    w("    type = PorousFlowPointSourceFromPostprocessor")
                    w("    variable = %s" % var)
                    w("    mass_flux = q_%s_%s_c%d" % (nm, ph, ci))
                    w("    point = '%.3f %.3f %.3f'" % tuple(c["xyz"]))
                    w("  []")
        w("[]")
        w("")
    w("[Postprocessors]")
    for aq in p["aquifers"]:
        w("  [aq%d_rate_kgps]   # + = water leaving the reservoir into aquifer %d" % (aq["id"], aq["id"]))
        w("    type = NodalSum")
        w("    variable = aq%d_flux" % aq["id"])
        w("    execute_on = 'initial timestep_end'")
        w("  []")
        w("  [aq%d_cum_out_kg]" % aq["id"])
        w("    type = TimeIntegratedPostprocessor")
        w("    value = aq%d_rate_kgps" % aq["id"])
        w("    time_integration_scheme = implicit-euler")
        w("    execute_on = 'initial timestep_end'")
        w("  []")
        w("  [aq%d_pressure_at_datum]" % aq["id"])
        w("    type = FunctionValuePostprocessor")
        w("    function = p_aq%d_fn" % aq["id"])
        w("    point = '0 0 %s'" % fmt(-aq["datum_depth"]))
        w("    execute_on = 'initial timestep_end'")
        w("  []")
    w("  [gas_in_place_kg]")
    w("    type = PorousFlowFluidMass")
    w("    fluid_component = 1")
    w("  []")
    w("  [water_in_place_kg]")
    w("    type = PorousFlowFluidMass")
    w("    fluid_component = 0")
    w("  []")
    for wl in p["wells"]:
        nm = safe(wl["name"])
        for ph, q in (("gas", wl["qg"]), ("water", wl["qw"])):
            if not any(q):
                continue
            w("  [%s_%s_rate_kgps]" % (nm, ph))
            w("    type = FunctionValuePostprocessor")
            w("    function = %s_%s_rate" % (nm, ph))
            w("    execute_on = 'initial timestep_begin'")
            w("  []")
            for ci, c in enumerate(wl["conns"]):
                w("  [q_%s_%s_c%d]" % (nm, ph, ci))
                w("    type = FunctionValuePostprocessor")
                w("    function = %s_%s_c%d" % (nm, ph, ci))
                w("    execute_on = 'initial timestep_begin'")
                w("    outputs = none")
                w("  []")
        c0 = wl["conns"][0]
        w("  [%s_pressure]" % nm)
        w("    type = PointValue")
        w("    variable = pwater")
        w("    point = '%.3f %.3f %.3f'" % tuple(c0["xyz"]))
        w("  []")
    w("[]")
    w("")
    w("[Dampers]   # limit Newton updates: robust phase appearance / disappearance")
    w("  [limit_p]")
    w("    type = MaxIncrement")
    w("    variable = pwater")
    w("    max_increment = 1e6")
    w("  []")
    w("  [limit_s]")
    w("    type = MaxIncrement")
    w("    variable = sgas")
    w("    max_increment = 0.2")
    w("  []")
    w("[]")
    w("")
    w("[Preconditioning]")
    w("  active = mumps")
    w("  [asm]")
    w("    type = SMP")
    w("    full = true")
    w("    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'")
    w("    petsc_options_value = ' asm      ilu          NONZERO                   2'")
    w("  []")
    w("  [mumps]   # direct solver, robust for small/medium models")
    w("    type = SMP")
    w("    full = true")
    w("    petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'")
    w("    petsc_options_value = ' lu       mumps'")
    w("  []")
    w("[]")
    w("")
    w("[Executioner]")
    w("  type = Transient")
    w("  solve_type = NEWTON")
    w("  line_search = none")
    w("  automatic_scaling = true")
    w("  end_time = %s # s (%.1f days)" % (fmt(p["end_time"] if p["end_time"] > 0 else 86400.0 * 365),
                                           (p["end_time"] if p["end_time"] > 0 else 86400.0 * 365) / 86400.0))
    w("  dtmax = %s # 30 days" % fmt(30 * 86400.0))
    if cap:
        w("  # with capillary pressure the model sits near capillary-gravity equilibrium: upwind directions flip on")
        w("  # links with almost no flow, and the residual stalls at |R| ~ 1-5 whatever its starting value")
        w("  nl_rel_tol = 1e-4")
        w("  nl_abs_tol = 10")
    else:
        w("  nl_rel_tol = 1e-6")
        w("  nl_abs_tol = 1e-3")
    w("  nl_max_its = 20")
    w("  l_max_its = 200")
    w("  [TimeStepper]")
    w("    type = IterationAdaptiveDT")
    w("    dt = 3600")
    w("    growth_factor = 1.5")
    w("    cutback_factor = 0.5")
    w("    optimal_iterations = 8")
    w("  []")
    w("[]")
    w("")
    w("[Outputs]")
    w("  [exodus]")
    w("    type = Exodus")
    w("  []")
    w("  [csv]")
    w("    type = CSV")
    times = rate_sync_times(p["wells"], p["end_time"])
    if len(times) > 3000:
        times = times[:: int(math.ceil(len(times) / 3000.0))]
    if times:
        w("    sync_times = '%s'" % wrap_list(times, 6, "                  "))
    w("  []")
    w("  print_linear_residuals = false")
    w("[]")
    return "\n".join(L) + "\n"


def moose_input_fv(p, args):
    """Cell-centred finite-volume PorousFlow input (two-point fluxes, like Eclipse)."""
    L = []
    w = L.append
    eq = p["eq"]

    def blk(name, typ, params, comment=""):
        w("  [%s]%s" % (name, ("   # " + comment) if comment else ""))
        w("    type = %s" % typ)
        for k_, v_ in params:
            w("    %s = %s" % (k_, v_))
        w("  []")

    w("# =============================================================================================")
    w("# MOOSE PorousFlow model generated by petrel2moose.py  (cell-centred finite volume, AD)")
    w("#   2-phase immiscible: phase 0 = water, phase 1 = gas; isothermal")
    w("#   primary variables (one value per cell): pwater [Pa], sgas [-]")
    w("#   mesh: %s   (x,y shifted by origin %.1f %.1f;  z = elevation, m)" % (p["mesh"], p["origin"][0], p["origin"][1]))
    w("#   Run:  porous_flow-opt -i <this file>   (or mpiexec -n N porous_flow-opt -i <this file>)")
    w("# =============================================================================================")
    w("")
    w("[Mesh]")
    blk("file", "FileMeshGenerator", [("file", p["mesh"]), ("use_for_exodus_restart", "true")])
    if p.get("nnc_file"):
        w("  allow_renumbering = false # the non-neighbour connections refer to the element ids of the file")
    w("[]")
    w("")
    w("[GlobalParams]")
    w("  PorousFlowDictator = dictator")
    w("  gravity = '0 0 -9.81'")
    w("[]")
    w("")
    w("[UserObjects]")
    blk("dictator", "PorousFlowDictator", [("porous_flow_vars", "'pwater sgas'"), ("number_fluid_phases", "2"),
                                          ("number_fluid_components", "2")])
    cap = p.get("capillary") == "deck" and any("pc" in rp for rp in p["relperm"])
    if not cap:
        blk("pc", "PorousFlowCapillaryPressureConst", [("pc", "0")])
    else:
        for rp in p["relperm"]:
            if "pc" in rp:
                c = rp["pc"]
                blk("pc_%s" % rp["block"], "PorousFlowCapillaryPressureBC",
                    [("pe", fmt(round(c["pe"], 2))), ("lambda", fmt(round(c["lam"], 4))),
                     ("sat_lr", fmt(round(c["sat_lr"], 4))), ("pc_max", fmt(round(c["pc_max"], 1)))],
                    "Brooks-Corey fit to SATNUM table %d, Pc(Swl) %s Pa" % (rp["table"], fmt(round(c["pcw"], 1))))
            else:
                blk("pc_%s" % rp["block"], "PorousFlowCapillaryPressureConst", [("pc", "0")])
    w("[]")
    w("")
    w("[Variables]")
    for v, f in (("pwater", "p_init"), ("sgas", "sg_init")):
        w("  [%s]" % v)
        w("    type = MooseVariableFVReal")
        if v == "pwater" and p["aquifers"]:
            w("    two_term_boundary_expansion = false # aquifer flux uses the cell pressure, as Eclipse")
        if p["init_mode"] == "restart":
            w("    initial_from_file_var = %s" % f)
        w("  []")
    w("[]")
    w("")

    fn, ic = [], []
    if p["init_mode"] != "restart":
        gwc_z, pg = 0.5 * (p["zmin"] + p["zmax"]), 1.5e7
        fn.append("  # !!! No EQUIL / restart available: EDIT gwc_elev and p_gwc below before running !!!")
        fn += ["  [p_ic]", "    type = ParsedFunction", "    symbol_names = 'gwc_elev p_gwc rho_w rho_g g'",
               "    symbol_values = '%s %s %s %s 9.81'" % (fmt(gwc_z), fmt(pg), fmt(p["rho_w"]), fmt(p["rho_g_res"] or 150.0)),
               "    expression = 'if(z > gwc_elev, p_gwc - rho_g * g * (z - gwc_elev), p_gwc - rho_w * g * (z - gwc_elev))'",
               "  []", "  [sg_ic]", "    type = ParsedFunction", "    symbol_names = 'gwc_elev sg_max'",
               "    symbol_values = '%s %s'" % (fmt(gwc_z), fmt(1.0 - p["swr"] - 1e-3)),
               "    expression = 'if(z > gwc_elev, sg_max, 0)'", "  []"]
        ic += ["  [p]", "    type = FunctionIC", "    variable = pwater", "    function = p_ic", "  []",
               "  [sg]", "    type = FunctionIC", "    variable = sgas", "    function = sg_ic", "  []"]
    for aq in p["aquifers"]:
        a = aq["id"]
        if aq.get("capacity"):
            fn += ["  [p_aq%d_fn]   # finite (Fetkovich) aquifer %d: hydrostatic about datum, minus cumulative "
                   "influx / (rho V ct)" % (a, a),
                   "    type = ParsedFunction", "    symbol_names = 'p0 datum_elev rho_w g cum cap'",
                   "    symbol_values = '%s %s %s 9.81 aq%d_cum_in_kg %s'"
                   % (fmt(aq["p0"]), fmt(-aq["datum_depth"]), fmt(p["rho_w"]), a, fmt(aq["capacity"])),
                   "    expression = 'p0 - rho_w * g * (z - datum_elev) - cum / cap'", "  []"]
        else:
            fn += ["  [p_aq%d_fn]   # aquifer %d pressure, hydrostatic about its datum (constant in time)" % (a, a),
                   "    type = ParsedFunction", "    symbol_names = 'p0 datum_elev rho_w g'",
                   "    symbol_values = '%s %s %s 9.81'" % (fmt(aq["p0"]), fmt(-aq["datum_depth"]), fmt(p["rho_w"])),
                   "    expression = 'p0 - rho_w * g * (z - datum_elev)'", "  []"]
    for wl in p["wells"]:
        nm = safe(wl["name"])
        for ph, q in (("gas", wl["qg"]), ("water", wl["qw"])):
            if not any(q):
                continue
            fn += ["  [%s_%s_rate]   # kg/s, + injection / - production" % (nm, ph),
                   "    type = PiecewiseConstant", "    direction = LEFT",
                   "    x = '%s'" % wrap_list(wl["t"], 6, "         "),
                   "    y = '%s'" % wrap_list(q, 6, "         "), "  []"]
    if fn:
        w("[Functions]")
        L.extend(fn)
        w("[]")
        w("")
    if ic:
        w("[ICs]")
        L.extend(ic)
        w("[]")
        w("")

    w("[AuxVariables]")
    for nm in ("poro", "perm_xx", "perm_yy", "perm_zz", "perm_xy"):
        w("  [%s]" % nm)
        w("    type = MooseVariableFVReal")
        w("    initial_from_file_var = %s" % nm)
        w("  []")
    for wl in p["wells"]:
        w("  [%s]   # well %s: rate fraction / cell volume in its completion cells [1/m3]" % (wl["weight_var"], wl["name"]))
        w("    type = MooseVariableFVReal")
        w("    initial_from_file_var = %s" % wl["weight_var"])
        w("  []")
    cwells = conn_wells(p, args)
    wwater = bool(cwells) and getattr(args, "well_water", "on") == "on"
    for wl in cwells:
        w("  [%s]   # well %s: connection factor / cell volume in its completion cells [1/m2] (COMPDAT CF in SI)" % (wl["cf_var"], wl["name"]))
        w("    type = MooseVariableFVReal")
        w("    initial_from_file_var = %s" % wl["cf_var"])
        w("  []")
        w("  [pi_%s]   # CF rho krg / mu_g of the completion cells, from the last time step [kg/s/Pa/m3]" % safe(wl["name"]))
        w("    type = MooseVariableFVReal")
        w("  []")
        if wwater:
            w("  [piw_%s]   # CF rho_w krw / mu_w of the completion cells, from the last time step [kg/s/Pa/m3]" % safe(wl["name"]))
            w("    type = MooseVariableFVReal")
            w("  []")
    if cwells:
        w("  [sg_begin]   # gas saturation at the start of the time step (time-step control at the wells)")
        w("    type = MooseVariableFVReal")
        w("  []")
        for nm, cm in (("krg_lag", "gas rel-perm"), ("mug_lag", "gas viscosity"), ("rho_g_lag", "gas density"))\
                + ((("krw_lag", "water rel-perm"), ("muw_lag", "water viscosity"), ("rho_w_lag", "water density"))
                   if wwater else ()):
            w("  [%s]   # %s at the start of the time step (well connections)" % (nm, cm))
            w("    type = MooseVariableFVReal")
            w("  []")
    for nm, val, cm in (("mass_frac_phase0_species0", 1, "water phase is pure water"),
                        ("mass_frac_phase1_species0", 0, "gas phase is pure gas")):
        w("  [%s]" % nm)
        w("    type = MooseVariableFVReal")
        w("    initial_condition = %d # %s" % (val, cm))
        w("  []")
    for aq in p["aquifers"]:
        w("  [p_aq%d]   # aquifer %d far-field pressure [Pa]" % (aq["id"], aq["id"]))
        w("    type = MooseVariableFVReal")
        w("    two_term_boundary_expansion = false # aquifer pressure at the cell depth, as Eclipse")
        w("  []")
    for nm in ("swater", "gas_density", "pgas"):
        w("  [%s]" % nm)
        w("    type = MooseVariableFVReal")
        w("  []")
    nnc = p.get("nnc_file")
    if nnc:
        for nm, cm in (("nnc_kr0", "water rel-perm"), ("nnc_kr1", "gas rel-perm"), ("nnc_mu0", "water viscosity"),
                       ("nnc_mu1", "gas viscosity"), ("nnc_rho0", "water density"), ("nnc_rho1", "gas density"),
                       ("nnc_pg", "gas pressure"), ("nnc_pc", "gas minus water pressure")):
            w("  [%s]   # %s at the start of the time step (non-neighbour connections)" % (nm, cm))
            w("    type = MooseVariableFVReal")
            w("  []")
    hyst = any("imb" in rp.get("tab", {}) for rp in p["relperm"])
    if hyst:
        w("  [sg_max]   # largest gas saturation each cell has had (rel-perm hysteresis); starts at the initial Sg")
        w("    type = MooseVariableFVReal")
        w("    initial_from_file_var = sg_init")
        w("  []")
    w("[]")
    w("")
    w("[AuxKernels]")
    if hyst:
        blk("sg_max", "ParsedAux", [("variable", "sg_max"), ("coupled_variables", "'sg_max sgas'"),
                                    ("expression", "'max(sg_max, sgas)'"), ("execute_on", "timestep_end")],
            "updated once per time step (lagged), as MOOSE's own hysteresis materials do")
    for aq in p["aquifers"]:
        blk("p_aq%d" % aq["id"], "FunctionAux", [("variable", "p_aq%d" % aq["id"]), ("function", "p_aq%d_fn" % aq["id"]),
                                                  ("execute_on", "'initial timestep_begin'")])
    blk("swater", "ADPorousFlowPropertyAux", [("variable", "swater"), ("property", "saturation"), ("phase", "0"),
                                               ("execute_on", "'initial timestep_end'")])
    blk("gas_density", "ADPorousFlowPropertyAux", [("variable", "gas_density"), ("property", "density"), ("phase", "1"),
                                                    ("execute_on", "'initial timestep_end'")])
    blk("pgas", "ADPorousFlowPropertyAux", [("variable", "pgas"), ("property", "pressure"), ("phase", "1"),
                                             ("execute_on", "'initial timestep_end'")])
    if cwells:
        blk("sg_begin", "ParsedAux", [("variable", "sg_begin"), ("coupled_variables", "sgas"), ("expression", "sgas"),
                                       ("execute_on", "'initial timestep_begin'")])
        for nm, prop, ph in (("krg_lag", "relperm", 1), ("mug_lag", "viscosity", 1), ("rho_g_lag", "density", 1))\
                + ((("krw_lag", "relperm", 0), ("muw_lag", "viscosity", 0), ("rho_w_lag", "density", 0)) if wwater else ()):
            blk(nm, "ADPorousFlowPropertyAux", [("variable", nm), ("property", prop), ("phase", ph),
                                                ("execute_on", "'initial timestep_begin'")], "lagged: fixed during the step")
        for wl in cwells:
            if wwater:
                blk("piw_%s" % safe(wl["name"]), "ParsedAux",
                    [("variable", "piw_%s" % safe(wl["name"])), ("coupled_variables", "'%s rho_w_lag krw_lag muw_lag'" % wl["cf_var"]),
                     ("expression", "'%s * rho_w_lag * krw_lag / muw_lag'" % wl["cf_var"]), ("execute_on", "'initial timestep_begin'")])
            blk("pi_%s" % safe(wl["name"]), "ParsedAux",
                [("variable", "pi_%s" % safe(wl["name"])), ("coupled_variables", "'%s rho_g_lag krg_lag mug_lag'" % wl["cf_var"]),
                 ("expression", "'%s * rho_g_lag * krg_lag / mug_lag'" % wl["cf_var"]), ("execute_on", "'initial timestep_begin'")])
    if nnc:
        for nm, prop, ph in (("nnc_kr0", "relperm", 0), ("nnc_kr1", "relperm", 1), ("nnc_mu0", "viscosity", 0),
                             ("nnc_mu1", "viscosity", 1), ("nnc_rho0", "density", 0), ("nnc_rho1", "density", 1),
                             ("nnc_pg", "pressure", 1)):
            blk(nm, "ADPorousFlowPropertyAux", [("variable", nm), ("property", prop), ("phase", ph),
                                                ("execute_on", "'initial timestep_begin'")])
        blk("nnc_pc", "ParsedAux", [("variable", "nnc_pc"), ("coupled_variables", "'nnc_pg pwater'"),
                                    ("expression", "'nnc_pg - pwater'"), ("execute_on", "'initial timestep_begin'")])
    w("[]")
    w("")
    w("[FVKernels]")
    if nnc:
        for comp, var, nm, xf in ((0, "pwater", "water", "'1 0'"), (1, "sgas", "gas", "'0 1'")):
            blk("nnc_%s" % nm, "FVPorousFlowNNCFlux",
                [("variable", var), ("nnc_file", nnc), ("pressure", "pwater"), ("pressure_offsets", "'0 nnc_pc'"),
                 ("relperm", "'nnc_kr0 nnc_kr1'"), ("density", "'nnc_rho0 nnc_rho1'"),
                 ("viscosity", "'nnc_mu0 nnc_mu1'"), ("mass_fraction", xf)],
                "%s flow through the non-neighbour (fault) connections: implicit pressure, mobility from the "
                "start of the step" % nm)
    for comp, var, nm in ((0, "pwater", "water"), (1, "sgas", "gas")):
        blk("mass_%s_dot" % nm, "FVPorousFlowMassTimeDerivative", [("fluid_component", comp), ("variable", var)])
        blk("flux_%s" % nm, "FVPorousFlowAdvectiveFluxTPFA" if args.fv_flux == "tpfa" else "FVPorousFlowAdvectiveFlux",
            [("fluid_component", comp), ("variable", var)])
    for wl in p["wells"]:
        nm = safe(wl["name"])
        for ph, var, q in (("gas", "sgas", wl["qg"]), ("water", "pwater", wl["qw"])):
            if not any(q):
                continue
            blk("%s_%s" % (nm, ph), "FVCoupledForce", [("variable", var), ("v", "src_%s_%s" % (nm, ph))],
                "well %s, %s rate distributed over its completion cells" % (wl["name"], ph))
        if wwater and wl in cwells:
            blk("%s_water_prod" % nm, "FVCoupledForce", [("variable", "pwater"), ("v", "src_%s_water_prod" % nm)],
                "well %s: water produced with the gas, by water mobility at the same drawdown" % wl["name"])
    w("[]")
    w("")
    if p["aquifers"]:
        w("[FVBCs]")
        for aq in p["aquifers"]:
            blk("aquifer_%d" % aq["id"], "FVFunctorNeumannBC",
                [("variable", "pwater"), ("boundary", "'%s'" % " ".join(aq["boundaries"])),
                 ("functor", "aq%d_influx" % aq["id"])],
                "aquifer %d (%s): water influx = C (p_aq - pwater)" % (aq["id"], aq["type"]))
        w("[]")
        w("")
    w("[FluidProperties]")
    if p.get("pvtw"):
        pw = p["pvtw"]
        rho_ref = p["rho_w_sc"] / pw["bw"]
        blk("water", "SimpleFluidProperties", [("density0", fmt(rho_ref * math.exp(-pw["cw"] * pw["pref"]))),
                                                ("bulk_modulus", fmt(1.0 / pw["cw"])), ("viscosity", fmt(pw["mu"])),
                                                ("thermal_expansion", "0")],
            "from PVTW: rho = rho_sc/Bw exp(cw (p - pref)), constant viscosity")
    else:
        blk("water", "Water97FluidProperties", [])
    if p.get("gas_pvt_file"):
        blk("gas", "TabulatedBicubicFluidProperties", [("fluid_property_file", p["gas_pvt_file"]),
                                                       ("interpolated_properties", "'density viscosity'"),
                                                       ("out_of_bounds_behavior", "set_to_closest_bound")],
            "from PVDG (Eclipse dry-gas table)")
    else:
        blk("gas", "MethaneFluidProperties", [])
    w("[]")
    w("")
    w("[Materials]")
    blk("temperature", "ADPorousFlowTemperature", [("temperature", "%s # K" % fmt(p["T"]))])
    if not cap:
        blk("ppss", "ADPorousFlow2PhasePS", [("phase0_porepressure", "pwater"), ("phase1_saturation", "sgas"),
                                             ("capillary_pressure", "pc")])
    else:
        for rp in p["relperm"]:
            blk("ppss_%s" % rp["block"], "ADPorousFlow2PhasePS",
                [("block", rp["block"]), ("phase0_porepressure", "pwater"), ("phase1_saturation", "sgas"),
                 ("capillary_pressure", "pc_%s" % rp["block"])])
    blk("massfrac", "ADPorousFlowMassFraction",
        [("mass_fraction_vars", "'mass_frac_phase0_species0 mass_frac_phase1_species0'")])
    for ph, nm in ((0, "water"), (1, "gas")):
        blk("%s" % nm, "ADPorousFlowSingleComponentFluid", [("fp", nm), ("phase", ph), ("compute_enthalpy", "false"),
                                                            ("compute_internal_energy", "false")])
    multi = len(p["relperm"]) > 1
    for rp in p["relperm"]:
        if rp.get("tab"):
            _write_relperm_tables(blk, rp, multi)
            continue
        for ph, nm in ((0, "water"), (1, "gas")):
            prm = []
            if multi:
                prm.append(("block", rp["block"]))
            sw_res = rp.get("swr_water", rp["swr"])
            prm += [("phase", ph), ("n", fmt(round(rp["nw"] if ph == 0 else rp["ng"], 3))),
                    ("s_res", fmt(round(sw_res if ph == 0 else rp["sgr"], 4))),
                    ("sum_s_res", fmt(round((sw_res if ph == 0 else rp.get("swr_gas", rp["swr"])) + rp["sgr"], 4))),
                    ("scaling", "%s # end point (SATNUM table %d)" % (fmt(round(rp["kw_max"] if ph == 0 else rp["kg_max"], 4)), rp["table"]))]
            blk("relperm_%s%s" % (nm, "_" + rp["block"] if multi else ""), "ADPorousFlowRelativePermeabilityCorey", prm)
    if p.get("rock"):
        # Eclipse ROCK: porosity = phi_ref (1 + X + X^2/2), X = c_r (P - p_ref), P the effective fluid pressure
        blk("eff_fluid_pressure", "ADPorousFlowEffectiveFluidPressure", [])
        blk("porosity", "ADPorousFlowPorosityRock",
            [("porosity", "poro"), ("reference_pressure", fmt(round(p["rock"]["p_ref"], 1))),
             ("compressibility", fmt(p["rock"]["c_r"]))],
            "Eclipse ROCK compressibility (add-on, moose_addon_fv); deck porosity at the reference pressure")
    else:
        blk("porosity", "ADPorousFlowPorosityConst", [("porosity", "poro")])
    blk("permeability", "ADPorousFlowPermeabilityConstFromVar",
        [("perm_xx", "perm_xx"), ("perm_yy", "perm_yy"), ("perm_zz", "perm_zz"), ("perm_xy", "perm_xy"), ("perm_yx", "perm_xy")])
    w("[]")
    w("")
    w("[FunctorMaterials]")
    # functor materials: well sources, aquifer influx, gas mass density for in-place reporting
    for wl in p["wells"]:
        nm = safe(wl["name"])
        for ph, q in (("gas", wl["qg"]), ("water", wl["qw"])):
            if not any(q):
                continue
            if ph == "gas" and wl in cwells:
                blk("cap_%s" % nm, "ADParsedFunctorMaterial",
                    [("property_name", "cap_%s" % nm), ("functor_names", "'pi_%s pwater rho_g_lag'" % nm),
                     ("functor_symbols", "'a p rho'"),
                     ("expression", "'a * max(p - 101325 - rho * 9.81 * (%s - z), 0)'" % fmt(wl["z_ref"]))],
                    "kg/s/m3 the cell can deliver with BHP = 1 atm at the reference depth %.2f m" % -wl["z_ref"])
                blk("src_%s_gas" % nm, "ADParsedFunctorMaterial",
                    [("property_name", "src_%s_gas" % nm),
                     ("functor_names", "'%s %s_gas_rate pi_%s sumpi_%s cap_%s'" % (wl["weight_var"], nm, nm, nm, nm)),
                     ("functor_symbols", "'w q a sa cap'"),
                     ("expression", "'if(q < 0, -min(-q * if(sa > 0, a / sa, w), cap), w * q)'")],
                    "kg/s/m3: production split by CF x mobility and limited by the cell's deliverability")
                if wwater:
                    blk("src_%s_water_prod" % nm, "ADParsedFunctorMaterial",
                        [("property_name", "src_%s_water_prod" % nm),
                         ("functor_names", "'%s_gas_rate sumpi_%s piw_%s pwater rho_g_lag'" % (nm, nm, nm)),
                         ("functor_symbols", "'q sa b p rho'"),
                         ("expression", "'if(q < 0 & sa > 0, -b * min(-q / sa, max(p - 101325 - rho * 9.81 * (%s - z), 0)), 0)'"
                          % fmt(wl["z_ref"]))],
                        "kg/s/m3 of water: the completion's drawdown for the gas (target / sum of gas mobilities, at "
                        "most down to BHP = 1 atm) times its water mobility CF rho_w krw / mu_w")
                continue
            blk("src_%s_%s" % (nm, ph), "ADParsedFunctorMaterial",
                [("property_name", "src_%s_%s" % (nm, ph)), ("functor_names", "'%s %s_%s_rate'" % (wl["weight_var"], nm, ph)),
                 ("functor_symbols", "'w q'"), ("expression", "'w * q'")], "kg/s/m3")
    for aq in p["aquifers"]:
        if aq.get("capacity"):
            tau = aq["capacity"] / (aq["mult"] * aq["area"])
            blk("aq%d_influx" % aq["id"], "ADParsedFunctorMaterial",
                [("property_name", "aq%d_influx" % aq["id"]), ("functor_names", "'p_aq%d pwater dt_size'" % aq["id"]),
                 ("functor_symbols", "'paq p dt'"),
                 ("expression", "'%s * %s / max(dt, 1) * (1 - exp(-max(dt, 1) / %s)) * (paq - p)'"
                  % (fmt(aq["mult"]), fmt(tau), fmt(tau)))],
                "kg/s/m2 into the reservoir; Eclipse's Fetkovich step average (aquifer time constant %.3g s)" % tau)
        else:
            blk("aq%d_influx" % aq["id"], "ADParsedFunctorMaterial",
                [("property_name", "aq%d_influx" % aq["id"]), ("functor_names", "'p_aq%d pwater'" % aq["id"]),
                 ("functor_symbols", "'paq p'"), ("expression", "'%s * (paq - p)'" % fmt(aq["mult"]))],
                "kg/s/m2 into the reservoir")
    blk("gas_mass_density", "ADParsedFunctorMaterial",
        [("property_name", "gas_mass_density"), ("functor_names", "'poro sgas gas_density'"),
         ("functor_symbols", "'phi s rho'"), ("expression", "'phi * s * rho'")], "kg of gas per m3 of rock")
    if cwells:
        blk("dsg_wells", "ADParsedFunctorMaterial",
            [("property_name", "dsg_wells"), ("functor_names", "'%s sgas sg_begin'" % " ".join(wl["weight_var"] for wl in cwells)),
             ("functor_symbols", "'%s s s0'" % " ".join("c%d" % i for i in range(len(cwells)))),
             ("expression", "'if(%s > 0, abs(s - s0), 0)'" % " + ".join("c%d" % i for i in range(len(cwells))))],
            "change of gas saturation during the step in the completion cells (weights ww > 0 there; the CF field is "
            "too small for the parser's comparison tolerance)")
    for wl in cwells:
        nm = safe(wl["name"])
        blk("wbp_%s" % nm, "ADParsedFunctorMaterial",
            [("property_name", "wbp_%s" % nm), ("functor_names", "'%s pgas gas_density'" % wl["cf_var"]),
             ("functor_symbols", "'c pg rho'"), ("expression", "'c * (pg + rho * 9.81 * (z - (%s)))'" % fmt(wl["z_ref"]))],
            "CF-weighted gas pressure corrected to the BHP reference depth (for WBP)")
    w("[]")
    w("")
    w("[Postprocessors]")
    if any(aq.get("capacity") for aq in p["aquifers"]):
        blk("dt_size", "TimestepSize", [("execute_on", "'initial timestep_begin'")])
    blk("gas_in_place_kg", "ElementIntegralFunctorPostprocessor", [("functor", "gas_mass_density"),
                                                                   ("execute_on", "'initial timestep_end'")])
    for aq in p["aquifers"]:
        blk("aq%d_influx_kgps" % aq["id"], "SideIntegralFunctorPostprocessor",
            [("functor", "aq%d_influx" % aq["id"]), ("boundary", "'%s'" % " ".join(aq["boundaries"])),
             ("execute_on", "'initial timestep_end'")])
        blk("aq%d_cum_in_kg" % aq["id"], "TimeIntegratedPostprocessor",
            [("value", "aq%d_influx_kgps" % aq["id"]), ("time_integration_scheme", "implicit-euler"),
             ("execute_on", "'initial timestep_end'")])
        blk("aq%d_pressure_at_datum" % aq["id"], "FunctionValuePostprocessor",
            [("function", "p_aq%d_fn" % aq["id"]), ("point", "'0 0 %s'" % fmt(-aq["datum_depth"])),
             ("execute_on", "'initial timestep_end'")])
    if cwells:
        blk("dsg_well_max", "ElementExtremeFunctorValue", [("functor", "dsg_wells"), ("value_type", "max"),
                                                          ("execute_on", "'initial timestep_end'")])
        blk("dt_taken", "TimestepSize", [("execute_on", "'initial timestep_end'"), ("outputs", "none")])
        blk("dt_well_limit", "ParsedPostprocessor",
            [("expression", "'max(dt_taken, 60) * %g / max(dsg_well_max, %g)'" % (args.well_dsg, args.well_dsg / 4)),
             ("pp_names", "'dt_taken dsg_well_max'"), ("execute_on", "'initial timestep_end'")],
            "the well mobility is taken from the start of the step: keep the saturation change of the completion cells small")
    for wl in cwells:
        nm = safe(wl["name"])
        blk("sumpi_%s" % nm, "ElementIntegralVariablePostprocessor",
            [("variable", "pi_%s" % nm), ("execute_on", "'initial timestep_begin'"), ("force_postaux", "true")],
            "sum over the completions of CF rho krg / mu_g (allocation of the target rate)")
        blk("%s_gas_rate_actual_kgps" % nm, "ElementIntegralFunctorPostprocessor",
            [("functor", "src_%s_gas" % nm), ("execute_on", "'initial timestep_end'")], "+ injection / - production")
        blk("%s_cum_gas_actual_kg" % nm, "TimeIntegratedPostprocessor",
            [("value", "%s_gas_rate_actual_kgps" % nm), ("execute_on", "'initial timestep_end'")])
        if wwater:
            blk("%s_water_rate_actual_kgps" % nm, "ElementIntegralFunctorPostprocessor",
                [("functor", "src_%s_water_prod" % nm), ("execute_on", "'initial timestep_end'")], "- production")
            blk("%s_cum_water_actual_kg" % nm, "TimeIntegratedPostprocessor",
                [("value", "%s_water_rate_actual_kgps" % nm), ("execute_on", "'initial timestep_end'")])
        blk("%s_cf_sum" % nm, "ElementIntegralVariablePostprocessor", [("variable", wl["cf_var"]), ("execute_on", "initial")])
        blk("%s_wbp_num" % nm, "ElementIntegralFunctorPostprocessor",
            [("functor", "wbp_%s" % nm), ("execute_on", "'initial timestep_end'"), ("outputs", "none")])
        blk("%s_wbp" % nm, "ParsedPostprocessor",
            [("expression", "'%s_wbp_num / %s_cf_sum'" % (nm, nm)), ("pp_names", "'%s_wbp_num %s_cf_sum'" % (nm, nm)),
             ("execute_on", "'initial timestep_end'")],
            "Eclipse-style WBP: CF-weighted gas pressure of the completion cells at the BHP reference depth (Pa)")
    for wl in p["wells"]:
        nm = safe(wl["name"])
        for ph, q in (("gas", wl["qg"]), ("water", wl["qw"])):
            if any(q):
                blk("%s_%s_rate_kgps" % (nm, ph), "FunctionValuePostprocessor",
                    [("function", "%s_%s_rate" % (nm, ph)), ("execute_on", "'initial timestep_end'")])
        blk("%s_pressure" % nm, "PointValue", [("variable", "pwater"), ("point", "'%.3f %.3f %.3f'" % tuple(wl["conns"][0]["xyz"]))],
            "cell pressure at the first completion")
        blk("%s_gas_pressure" % nm, "PointValue", [("variable", "pgas"), ("point", "'%.3f %.3f %.3f'" % tuple(wl["conns"][0]["xyz"]))],
            "gas pressure in that cell (Eclipse's WBP is the gas pressure)")
    w("[]")
    w("")
    w("[Preconditioning]")
    w("  active = mumps")
    blk("mumps", "SMP", [("full", "true"), ("petsc_options_iname", "'-pc_type -pc_factor_mat_solver_type'"),
                         ("petsc_options_value", "' lu       mumps'")], "direct solver: robust, fine to ~1e6 unknowns")
    blk("asm", "SMP", [("full", "true"),
                       ("petsc_options_iname", "'-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap -ksp_gmres_restart'"),
                       ("petsc_options_value", "' asm      lu           NONZERO                   2               200'")],
        "parallel alternative for many cores")
    w("[]")
    w("")
    end = p["end_time"] if p["end_time"] > 0 else 86400.0 * 365
    w("[Executioner]")
    w("  type = Transient")
    w("  solve_type = NEWTON")
    w("  line_search = none")
    w("  automatic_scaling = true")
    w("  end_time = %s # s (%.1f days)" % (fmt(end), end / 86400.0))
    w("  dtmax = %s # 30 days" % fmt(30 * 86400.0))
    if cap:
        w("  nl_rel_tol = 1e-4 # near capillary-gravity equilibrium upwind directions flip on links with ~no flow")
        w("  nl_abs_tol = 10")
    else:
        w("  nl_rel_tol = 1e-6")
        w("  nl_abs_tol = 1e-3")
    w("  nl_max_its = 20")
    w("  l_max_its = 100")
    w("  [TimeStepper]")
    w("    type = IterationAdaptiveDT")
    w("    dt = 3600")
    w("    growth_factor = 1.5")
    w("    cutback_factor = 0.5")
    w("    optimal_iterations = 6")
    if conn_wells(p, args):
        w("    timestep_limiting_postprocessor = dt_well_limit # gas saturation of the completion cells changes <= %g per step" % args.well_dsg)
    w("  []")
    w("[]")
    w("")
    w("[Outputs]")
    w("  [exodus]")
    w("    type = Exodus")
    if getattr(args, "exodus_vars", "key") == "key":
        keyv = ["pwater", "pgas", "sgas", "swater"] + (["sg_max"] if hyst else [])
        w("    show = '%s' # what the analysis and the H2 generator need; --exodus-vars all writes every field" % " ".join(keyv))
    ed = getattr(args, "exodus_days", 10.0)
    if ed and ed > 0:
        w("    min_simulation_time_interval = %s # write at most every %g days (--exodus-days); the CSV has every step"
          % (fmt(ed * 86400.0), ed))
        w("    execute_on = 'initial timestep_end final'")
    w("  []")
    w("  [csv]")
    w("    type = CSV")
    times = rate_sync_times(p["wells"], end)
    if len(times) > 3000:
        times = times[:: int(math.ceil(len(times) / 3000.0))]
    if times:
        w("    sync_times = '%s'" % wrap_list(times, 6, "                  "))
    w("  []")
    w("  print_linear_residuals = false")
    w("[]")
    return "\n".join(L) + "\n"


def w_functions_wells(w, p):
    for wl in p["wells"]:
        nm = safe(wl["name"])
        for ph, q in (("gas", wl["qg"]), ("water", wl["qw"])):
            if not any(q):
                continue
            w("  [%s_%s_rate]   # kg/s, + injection / - production" % (nm, ph))
            w("    type = PiecewiseConstant")
            w("    direction = LEFT")
            w("    x = '%s'" % wrap_list(wl["t"], 6, "         "))
            w("    y = '%s'" % wrap_list(q, 6, "         "))
            w("  []")
            for ci, c in enumerate(wl["conns"]):
                w("  [%s_%s_c%d]" % (nm, ph, ci))
                w("    type = ParsedFunction")
                w("    symbol_names = 'q'")
                w("    symbol_values = '%s_%s_rate'" % (nm, ph))
                w("    expression = '%s * q'" % fmt(c["frac"]))
                w("  []")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_argument_group("input (one of)")
    src.add_argument("--case", help="Eclipse .DATA deck written by Petrel (EGRID/INIT/UNRST found next to it)")
    src.add_argument("--egrid", help=".EGRID file (if not next to the deck)")
    src.add_argument("--init", help=".INIT file with PORO/PERM*/NTG")
    src.add_argument("--unrst", help=".UNRST file for initial pressure/saturation")
    src.add_argument("--grdecl", help="ASCII GRDECL grid (Petrel: Export > Eclipse grid)")
    src.add_argument("--props", nargs="*", help="extra GRDECL files with PORO, PERMX, ... (global size)")
    ap.add_argument("-o", "--output", default="model", help="output prefix (default: model)")
    ap.add_argument("--restart-step", type=int, default=None,
                    help="initialise from UNRST report step index (0 = first, -1 = last). Default: EQUIL")
    ap.add_argument("--average-rates", type=float, default=0,
                    help="average the well-rate history over periods of N days (mass conserving; fewer time steps)")
    ap.add_argument("--merge-layers", type=int, default=0,
                    help="vertical upscaling: merge every N layers (recommended for thin-layered geomodels)")
    ap.add_argument("--mobile-connate-water", action="store_true",
                    help="keep table water residuals everywhere (water above them in the gas cap then drains, Pc = 0)")
    ap.add_argument("--swi-bin", type=float, default=0.025, help="initial-Sw bin width for the immobile-water sub-blocks")
    ap.add_argument("--capillary", choices=("none", "deck"), default="none",
                    help="capillary pressure: none (Pc = 0 with immobile connate water, default) or deck (Brooks-Corey "
                         "fits to the deck's Pc tables, scaled cell by cell like Eclipse's PCW/SWLPC; initial water "
                         "pressure from the restart)")
    ap.add_argument("--pc-bin-ratio", type=float, default=1.2,
                    help="with --capillary deck: cells whose PCW differs by less than this factor share a Pc curve")
    ap.add_argument("--pc-max-factor", type=float, default=10.0,
                    help="with --capillary deck: Pc is capped (smoothly) at this multiple of the cell's PCW")
    ap.add_argument("--pc-swl-bin", type=float, default=0.02,
                    help="with --capillary deck: SWLPC bin width for the Pc sub-blocks")
    ap.add_argument("--relperm", choices=("deck", "table", "corey"), default="deck",
                    help="with --fv: deck = the deck's SWFN/SGFN tables, with Eclipse's rel-perm hysteresis when the deck "
                         "has SATOPTS HYSTER (default; needs the add-on ADPorousFlowRelativePermeabilityHysteresisTable); "
                         "table = the drainage tables only; corey = Corey fits (as before). The nodal input always "
                         "uses Corey fits")
    ap.add_argument("--well-model", choices=("connection", "rate"), default="connection",
                    help="with --fv: connection = gas production split over the completions by CF x gas mobility and "
                         "limited by what each cell can deliver at BHP 1 atm (default; needs COMPDAT CF); rate = the "
                         "target rate taken from the completion cells by Kh share whatever their state (as before)")
    ap.add_argument("--exodus-days", type=float, default=10.0,
                    help="--fv: write the Exodus file at most every this many days (0 = every time step); the CSV "
                         "keeps every step. The full Croft grid writes ~8 MB per field per output")
    ap.add_argument("--exodus-vars", choices=["key", "all"], default="key",
                    help="--fv: 'key' writes pwater, pgas, sgas, swater (and sg_max with hysteresis) to the Exodus "
                         "file; 'all' also the static and diagnostic fields (permeability, lagged mobilities, ...)")
    ap.add_argument("--well-water", choices=["on", "off"], default="on",
                    help="connection wells: also produce water from each completion in proportion to its water "
                         "mobility at the gas drawdown (as Eclipse does); the history gas rate stays the target")
    ap.add_argument("--well-dsg", type=float, default=0.05,
                    help="connection wells: target change of gas saturation per time step in the completion cells "
                         "(time-step control, like Eclipse's DSMAX; default 0.05)")
    ap.add_argument("--no-rock-compressibility", action="store_true",
                    help="with --fv: ignore the deck's ROCK compressibility")
    ap.add_argument("--fv-flux", choices=("tpfa", "moose"), default="tpfa",
                    help="with --fv: tpfa = FVPorousFlowAdvectiveFluxTPFA (harmonic transmissibility, hydrostatically "
                         "consistent; needs the add-on kernel), moose = MOOSE's FVPorousFlowAdvectiveFlux")
    ap.add_argument("--fv", action="store_true",
                    help="write a cell-centred finite-volume input instead of the default nodal finite-element one "
                         "(experimental: MOOSE's FV PorousFlow flux is not hydrostatically consistent on dipping grids)")
    ap.add_argument("--restart-pressure", action="store_true",
                    help="use the restart PRESSURE as is (default with EQUIL: hydrostatic pressure + restart saturations)")
    ap.add_argument("--zones", help="blocks by layer ranges, e.g. 'Waarre_C:1-12,Waarre_A:13-30'")
    ap.add_argument("--block-by", help="blocks from an integer region property, e.g. FIPNUM or SATNUM")
    ap.add_argument("--min-poro", type=float, default=1e-4, help="drop cells with PORO*NTG <= this")
    ap.add_argument("--merge-tol", type=float, default=1e-3, help="node merge tolerance on a pillar [m]")
    ap.add_argument("--fault-snap", type=float, default=None,
                    help="merge nodes across faults whose throw is below min(this [m], --fault-snap-frac x cell "
                         "thickness), so faces offset by a small throw are connected (Eclipse connects them); "
                         "0 = off. Default: 0 with --nnc (the NNCs connect the faults), else 1.0")
    ap.add_argument("--nnc", choices=["on", "off"], default="on",
                    help="--fv: carry Eclipse's connections across faults (TRANX/Y/Z, NNC/TRANNNC from the "
                         "EGRID/INIT) as non-neighbour connections wherever the MOOSE cells share no face "
                         "(needs FVPorousFlowNNCFlux in the MOOSE build; default on when EGRID and INIT exist)")
    ap.add_argument("--nnc-mult", type=float, default=1.0, help="multiplier on the NNC transmissibilities")
    ap.add_argument("--fault-snap-frac", type=float, default=0.3, help="see --fault-snap")
    ap.add_argument("--kvkh", type=float, default=0.1, help="kv/kh if PERMZ is missing")
    ap.add_argument("--no-ntg-perm", action="store_true", help="do not multiply PERMX/PERMY by NTG")
    ap.add_argument("--no-shift", action="store_true", help="keep absolute map coordinates")
    ap.add_argument("--temperature", type=float, default=80.0, help="reservoir T [C] if no RTEMP in deck")
    ap.add_argument("--zfactor", type=float, default=0.9, help="gas Z used only for the hydrostatic IC")
    ap.add_argument("--mu-w", type=float, default=5e-4, help="water viscosity [Pa s] for aquifer estimate")
    ap.add_argument("--aquifer-mult", type=float, default=1e-9, help="fallback aquifer multiplier [kg/s/m2/Pa]")
    ap.add_argument("--end-time", type=float, default=365.0, help="end time [days] when no schedule")
    args = ap.parse_args(argv)
    od = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(od, exist_ok=True)
    convert(args)


if __name__ == "__main__":
    main()

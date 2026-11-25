#!/usr/bin/env python
"""Write a finite-volume MOOSE PorousFlow input for cyclic hydrogen storage, starting from the end of a
finite-volume history-match run made with em2ex2moose.py --fv (any grid: 19 or 92 layers).

    python gen_h2_fv.py --hist-input croft_full_fv.i --hist-exodus croft_full_fv_exodus.e -o croft_h2_fv

Reads
  * the history input (.i): mesh name, aquifer(s), well(s), water PVT, capillary pressure, rel-perm,
    rock compressibility and permeability blocks are copied from it, so the storage model has exactly
    the history model's rock and grid;
  * the history Exodus output: pressure and gas saturation of every cell at the last time step, and the
    cumulative aquifer influx (to re-base the aquifer pressure to its depleted state);
  * h2_mix_pvt.csv: density and viscosity of the native-gas + H2 mixture (build_mix_table.py).
Writes
  <out>.e   the history mesh with p_init / sg_init replaced by the end-of-history state
  <out>.i   the storage input
  <out>_report.md
Physics: 3 components (0 water, 1 native gas, 2 hydrogen) in 2 phases (water, gas), primary variables
pwater, sgas, zh2 (H2 mass fraction of the gas phase). Mixing of H2 into the native gas by
advection and the numerical dispersion of the upwind scheme. See README.md.
"""
import argparse, os, re, shutil, sys
import numpy as np
import netCDF4
from scipy.spatial import cKDTree

DAY = 86400.0
M_H2 = 2.016  # g/mol


def names(ds, v):
    return ["".join(c.decode() for c in row if c) for row in ds[v][:]]


def section(txt, name):
    m = re.search(r"^\[%s\]\n(.*?)^\[\]" % re.escape(name), txt, re.S | re.M)
    if not m:
        raise SystemExit("history input has no [%s] section" % name)
    return m.group(1)


def subblocks(sec, pattern):
    """sub-blocks '  [name] ... \n  []' of a section whose name matches the regex pattern"""
    out = []
    for m in re.finditer(r"^  \[([^\]]+)\]([^\n]*)\n(.*?)^  \[\]\n", sec, re.S | re.M):
        if re.fullmatch(pattern, m.group(1)):
            out.append((m.group(1), m.group(0).rstrip("\n")))
    return out


def param(block, key):
    m = re.search(r"^\s*%s\s*=\s*('[^']*'|\S+)" % key, block, re.M)
    return m.group(1).strip("'") if m else None


def centroids(ds):
    X = np.c_[ds["coordx"][:], ds["coordy"][:], ds["coordz"][:]]
    nb = ds.dimensions["num_el_blk"].size
    C = [ds["connect%d" % (b + 1)][:] - 1 for b in range(nb)]
    return np.vstack([X[c].mean(1) for c in C]), [len(c) for c in C]


def elem_field(ds, name, step):
    ev = names(ds, "name_elem_var")
    k = ev.index(name)
    nb = ds.dimensions["num_el_blk"].size
    out = []
    for b in range(nb):
        v = "vals_elem_var%deb%d" % (k + 1, b + 1)
        n = ds["connect%d" % (b + 1)].shape[0]
        out.append(ds[v][step] if v in ds.variables else np.full(n, np.nan))
    return np.concatenate(out)


def copy_with_elem_vars(src, dst, new):
    """Copy an Exodus (netCDF-3) file, adding element variables `new` (names), initialised to 0.
    All variables are defined first and filled afterwards (netCDF-3 rewrites the file on every re-definition)."""
    a = netCDF4.Dataset(src)
    a.set_auto_mask(False)
    b = netCDF4.Dataset(dst, "w", format=a.file_format)
    b.set_auto_mask(False)
    b.setncatts({k: a.getncattr(k) for k in a.ncattrs()})
    nev = a.dimensions["num_elem_var"].size if "num_elem_var" in a.dimensions else 0
    for k, d in a.dimensions.items():
        b.createDimension(k, None if d.isunlimited() else (nev + len(new) if k == "num_elem_var" else d.size))
    if "num_elem_var" not in a.dimensions:
        b.createDimension("num_elem_var", len(new))
    nblk = a.dimensions["num_el_blk"].size
    copy = [k for k in a.variables if k not in ("name_elem_var", "elem_var_tab")]
    for k in copy:
        v = a.variables[k]
        o = b.createVariable(k, v.dtype, v.dimensions)
        o.setncatts({x: v.getncattr(x) for x in v.ncattrs()})
    b.createVariable("name_elem_var", "S1", ("num_elem_var", "len_name"))
    b.createVariable("elem_var_tab", "i4", ("num_el_blk", "num_elem_var"))
    for j in range(len(new)):
        for e in range(nblk):
            b.createVariable("vals_elem_var%deb%d" % (nev + j + 1, e + 1), "f8", ("time_step", "num_el_in_blk%d" % (e + 1)))
    for k in copy:   # data
        b[k][:] = a[k][:]
    ln = a.dimensions["len_name"].size
    rows = [] if nev == 0 else [np.asarray(r, "S1") for r in a["name_elem_var"][:]]
    rows += [np.array([c.encode() for c in n_.ljust(ln, "\0")], dtype="S1") for n_ in new]
    b["name_elem_var"][:] = np.array(rows)
    t0 = a["elem_var_tab"][:] if "elem_var_tab" in a.variables else np.ones((nblk, nev), "i4")
    b["elem_var_tab"][:] = np.hstack([t0, np.ones((nblk, len(new)), "i4")])
    nt = max(a.dimensions["time_step"].size, 1)
    for j in range(len(new)):
        for e in range(nblk):
            b["vals_elem_var%deb%d" % (nev + j + 1, e + 1)][:] = np.zeros((nt, a.dimensions["num_el_in_blk%d" % (e + 1)].size))
    a.close()
    b.close()


def write_zh2_table(src, dst):
    """The mixture table with its composition column (native-gas fraction xng, in the 'temperature' slot) turned
    into the H2 fraction zh2 = 1 - xng, rows sorted as TabulatedFluidProperties needs (pressure, then zh2)."""
    with open(src) as f:
        head = f.readline()
    t = np.loadtxt(src, delimiter=",", skiprows=1)
    t[:, 1] = 1.0 - t[:, 1]
    t = t[np.lexsort((t[:, 1], t[:, 0]))]
    with open(dst, "w", newline="\n") as f:
        f.write(head)
        for r in t:
            f.write(",".join("%.10g" % v for v in r) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hist-input", required=True, help="history-match FV input (.i) written by em2ex2moose.py --fv")
    ap.add_argument("--hist-exodus", required=True, help="Exodus output of that run (uses its last time step)")
    ap.add_argument("--mesh", help="mesh of the history run (default: the file named in the history input)")
    ap.add_argument("--mix-table", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "h2_mix_pvt.csv"))
    ap.add_argument("-o", "--out", default="croft_h2_fv", help="output base name")
    ap.add_argument("--mass-t", type=float, default=1000.0, help="H2 injected per cycle and gas withdrawn per cycle [t]")
    ap.add_argument("--inj-days", type=float, default=182.625)
    ap.add_argument("--wd-days", type=float, default=91.3125)
    ap.add_argument("--cycles", type=int, default=10)
    ap.add_argument("--dtmax-days", type=float, default=7.0)
    ap.add_argument("--eps", type=float, default=1e-5, help="dissolved-gas mass fraction in water (numerical)")
    ap.add_argument("--native-molar-mass", type=float, default=18.92, help="native gas molar mass [g/mol]")
    ap.add_argument("--exodus-days", type=float, default=30.0,
                    help="also write the 3D fields every this many days (default 30; 0: only at the switch "
                         "times). Each output of the full grid is about 50 MB")
    ap.add_argument("--deck", help="Eclipse .DATA deck: well connection factors (COMPDAT CF), BHP reference depth "
                    "(WELSPECS), initial pressure (EQUIL + PVDG). Without it the Croft values are used.")
    ap.add_argument("--p-max-psia", type=float, help="injection limit on the bottom-hole pressure at the BHP reference "
                    "depth (default: the initial, pre-production reservoir pressure at that depth)")
    ap.add_argument("--p-min-psia", type=float, help="withdrawal limit (minimum bottom-hole pressure) at the reference "
                    "depth, also used for the deliverability (default: the lowest well-block pressure of the history run)")
    ap.add_argument("--rho-sc-h2", type=float, default=0.08527, help="H2 density at standard conditions [kg/sm3] (15 C, 1 atm)")
    ap.add_argument("--nl-rel-tol", type=float, default=1e-3)
    ap.add_argument("--nl-abs-tol", type=float, default=1e-4, help="kg/s")
    ap.add_argument("--mumps-relax", type=int, default=100, help="MUMPS ICNTL(14) workspace relaxation [%%]")
    ap.add_argument("--mumps-ooc", action="store_true", help="MUMPS out-of-core factorisation (less memory, slower)")
    a = ap.parse_args()

    txt = open(a.hist_input).read().replace("\r\n", "\n")
    hdir = os.path.dirname(os.path.abspath(a.hist_input))
    rep = ["# Hydrogen storage input %s.i\n" % a.out, "History input: `%s`, history output: `%s`\n" % (a.hist_input, a.hist_exodus)]

    # ------------------------------------------------------------------ history input
    if "FVPorousFlowAdvectiveFlux" not in txt:
        raise SystemExit("the history input is not a finite-volume input (em2ex2moose.py --fv)")
    mesh = a.mesh or os.path.join(hdir, param(section(txt, "Mesh"), "file"))
    uo = section(txt, "UserObjects")
    pcs = subblocks(uo, r"pc\w*")
    mats = section(txt, "Materials")
    keep = subblocks(mats, r"ppss\w*|relperm_\w+|eff_fluid_pressure|porosity|permeability")
    water = [b for n, b in subblocks(section(txt, "FluidProperties"), "water")][0]
    funcs = section(txt, "Functions")
    fvbcs = section(txt, "FVBCs")
    fmats = section(txt, "FunctorMaterials")
    pps = section(txt, "Postprocessors")
    auxv = section(txt, "AuxVariables")
    aqs = []
    for n, b in subblocks(funcs, r"p_aq(\d+)_fn"):
        k = re.search(r"\d+", n).group(0)
        p0, zd, rho, g, cumname, cap = param(b, "symbol_values").split()
        bnd = param([bb for nn, bb in subblocks(fvbcs, "aquifer_%s" % k)][0], "boundary")
        expr = param([bb for nn, bb in subblocks(fmats, "aq%s_influx" % k)][0], "expression")
        aqs.append(dict(k=k, p0=float(p0), zd=zd, rho=rho, g=g, cap=float(cap), bnd=bnd, expr=expr))
    if not aqs:
        raise SystemExit("no aquifer found in the history input")
    wells = [re.match(r"ww_(.+)", n).group(1) for n, b in subblocks(auxv, r"ww_.+")]
    pts = {}
    for w in wells:
        b = [bb for nn, bb in subblocks(pps, r"%s_gas_pressure" % re.escape(w))]
        pts[w] = param(b[0], "point")
    hist_end = float(param(section(txt, "Executioner"), "end_time"))
    nnc_src = None
    for n_, b_ in subblocks(section(txt, "FVKernels"), r"nnc_\w+"):
        if "FVPorousFlowNNCFlux" in b_:
            nnc_src = os.path.join(hdir, param(b_, "nnc_file"))
            break

    # ------------------------------------------------------------------ end-of-history state
    he = netCDF4.Dataset(a.hist_exodus)
    he.set_auto_mask(False)
    times = he["time_whole"][:]
    step = len(times) - 1
    if abs(times[step] - hist_end) > 1.0:
        print("WARNING: last time in %s is %.0f s, history end_time is %.0f s" % (a.hist_exodus, times[step], hist_end))
    rep.append("End-of-history state from time step %d, t = %.0f s (%.1f days); history end_time %.0f s.\n"
               % (step, times[step], times[step] / DAY, hist_end))
    glo = names(he, "name_glo_var")
    gv = he["vals_glo_var"][:]
    for aq in aqs:
        cum = float(gv[step, glo.index("aq%s_cum_in_kg" % aq["k"])])
        aq["cum_end"] = cum
        aq["p0_new"] = aq["p0"] - cum / aq["cap"]
        rep.append("Aquifer %s: cumulative influx at the end of history %.4g kg (%.3f MMstb); datum pressure "
                   "re-based from %.0f Pa to %.0f Pa (%.1f psi lower).\n"
                   % (aq["k"], cum, cum / 1020.34 / 0.158987 / 1e6, aq["p0"], aq["p0_new"], (aq["p0"] - aq["p0_new"]) / 6894.757))
    gip_hist = float(gv[step, glo.index("gas_in_place_kg")]) if "gas_in_place_kg" in glo else None
    pw_h = elem_field(he, "pwater", step)
    sg_h = elem_field(he, "sgas", step)
    ch, _ = centroids(he)
    # rel-perm hysteresis (history input made with --relperm deck): carry each cell's largest gas saturation
    hyst = any("max_saturation" in b for n, b in keep if n.startswith("relperm_"))
    sgmax_h = None
    if hyst:
        if "sg_max" in names(he, "name_elem_var"):
            sgmax_h = elem_field(he, "sg_max", step)
            rep.append("Rel-perm hysteresis: each cell's largest gas saturation so far (sg_max) is taken from the "
                       "history output and keeps evolving in the storage run.\n")
        else:
            sgmax_h = np.maximum(elem_field(he, "sgas", 0), sg_h)
            rep.append("Rel-perm hysteresis: the history output has no sg_max (history run without hysteresis), so "
                       "each cell's largest gas saturation is taken as max(initial, end-of-history) Sg. Gas "
                       "saturation fell almost everywhere during depletion, so this is close to the true maximum; "
                       "cells already invaded by water start on their Killough scanning curves.\n")
            print("NOTE: no sg_max in the history output; sg_max = max(initial, end-of-history) gas saturation")

    # ------------------------------------------------------------------ new mesh with that state
    out_e = a.out + ".e"
    if hyst:
        copy_with_elem_vars(mesh, out_e, ["sg_max"])
    else:
        shutil.copyfile(mesh, out_e)
    me = netCDF4.Dataset(out_e, "a")
    me.set_auto_mask(False)
    cm, nper = centroids(me)
    dist, idx = cKDTree(ch).query(cm)
    if dist.max() > 1e-3 * max(1.0, np.ptp(cm[:, 2])):
        raise SystemExit("history output and mesh do not match (largest centroid distance %.3g m)" % dist.max())
    ev = names(me, "name_elem_var")
    off = 0
    for b, n in enumerate(nper):
        fields = [("p_init", pw_h), ("sg_init", np.clip(sg_h, 0.0, 0.999))]
        if hyst:
            fields.append(("sg_max", np.clip(sgmax_h, 0.0, 1.0)))
        for name, val in fields:
            v = "vals_elem_var%deb%d" % (ev.index(name) + 1, b + 1)
            me[v][0, :] = val[idx[off:off + n]]
        off += n
    poro = elem_field(me, "poro", 0)
    pv = elem_field(me, "pore_volume", 0)
    me.close()
    sg0 = np.clip(sg_h[idx], 0, 0.999)
    rep.append("Mesh `%s` copied to `%s`; %d cells matched to the history output (largest centroid distance %.2g m); "
               "p_init and sg_init now hold the end-of-history pwater and sgas.\n" % (mesh, out_e, len(cm), dist.max()))
    rep.append("Initial state: %d gas cells (Sg > 1e-3), gas pore volume %.4g m3, pwater %.3f-%.3f MPa.\n"
               % ((sg0 > 1e-3).sum(), (pv * sg0).sum(), pw_h.min() / 1e6, pw_h.max() / 1e6))
    if gip_hist is not None:
        rep.append("Gas in place at the end of history (history run): %.6g kg.\n" % gip_hist)

    # ------------------------------------------------------------------ schedule
    EPS = a.eps
    q_inj = a.mass_t * 1e3 / (a.inj_days * DAY)
    q_wd = a.mass_t * 1e3 / (a.wd_days * DAY)
    P = (a.inj_days + a.wd_days) * DAY
    Ti = a.inj_days * DAY
    t_end = a.cycles * P
    sw = []
    for k in range(a.cycles):
        sw += [k * P + Ti, (k + 1) * P]
    sw_s = " ".join("%.1f" % s for s in sw)
    ex_t = set(round(s, 1) for s in sw)
    if a.exodus_days > 0:
        ex_t |= set(round(t, 1) for t in np.arange(a.exodus_days * DAY, a.cycles * P, a.exodus_days * DAY))
    ex_s = " ".join("%.1f" % t for t in sorted(ex_t))
    sched_y = " ".join(str(i) for i in range(len(sw) + 1))
    well = wells[0]
    if len(wells) > 1:
        print("NOTE: %d wells in the history input; injection and withdrawal use %s" % (len(wells), well))
    pt = pts[well]
    MNG = a.native_molar_mass
    # ---------------------------------------------------------------- well model and pressure limits
    PSI_ = 6894.757293168
    CF_SI = 1e-3 * 0.158987294928 / 86400.0 / PSI_  # COMPDAT CF [cP rb/d/psi] -> [m3]
    if a.deck:
        _here = os.path.dirname(os.path.abspath(__file__))
        for _p in (os.environ.get("EM2EX_TO_MOOSE", ""), os.path.join(_here, "..", "em2ex_to_moose"), _here):
            if _p and os.path.exists(os.path.join(_p, "eclio.py")):
                sys.path.insert(0, _p)
                break
        import eclio
        dk = eclio.Deck(a.deck)
        sch = dk.schedule()
        cf_tot = sum(c["cf"] for c in sch["compdat"][well] if str(c.get("status", "OPEN")).upper() == "OPEN")
        z_ref_m = sch["welspecs"][well]["ref_depth"] * 0.3048
        eq = dk.kw["EQUIL"][0][0]
        pg_, bg_, _mu = dk.pvdg()
        rho_sc_ng = dk.densities()[2]
        p_dat = eq[1] * PSI_
        rho_ref = rho_sc_ng / np.interp(p_dat, pg_, bg_)
        p_init_ref = p_dat - rho_ref * 9.81 * (eq[0] * 0.3048 - z_ref_m)
        src_txt = "deck %s" % a.deck
    else:
        cf_tot, z_ref_m, rho_ref, rho_sc_ng = 3.473876, 6509.19 * 0.3048, 145.688, 0.7987
        p_init_ref = 2865.882 * PSI_
        src_txt = "Croft defaults (no --deck given)"
    C_well = cf_tot * CF_SI
    p_max = a.p_max_psia * PSI_ if a.p_max_psia else p_init_ref
    if a.p_min_psia:
        p_min = a.p_min_psia * PSI_
    else:
        if "%s_gas_pressure" % well not in glo:
            raise SystemExit("the history output has no %s_gas_pressure: give the minimum BHP with --p-min-psia" % well)
        k = glo.index("%s_gas_pressure" % well)
        hist_p = gv[1:, k]
        p_min = float(hist_p[hist_p > 0].min())
    # H2 density at reservoir conditions from the mixture table (xng = 0): rho = p / (A + B p)
    mt = np.loadtxt(a.mix_table, delimiter=",", skiprows=1)
    h2 = mt[np.isclose(mt[:, 1], 0.0)]
    sel = (h2[:, 0] > 5e6) & (h2[:, 0] < 35e6)
    Bh, Ah = np.polyfit(h2[sel, 0], h2[sel, 0] / h2[sel, 2], 1)
    rep.append("Well %s: sum of open COMPDAT CF = %.6g cP.rb/d/psi = %.6g m3 (%s); BHP reference depth %.2f m.\n"
               % (well, cf_tot, C_well, src_txt, z_ref_m))
    rep.append("Injection limit: BHP <= %.1f psia at the reference depth (initial reservoir pressure there: %.1f psia; "
               "gradient %.1f kg/m3). Withdrawal limit: BHP >= %.1f psia.\n"
               % (p_max / PSI_, p_init_ref / PSI_, rho_ref, p_min / PSI_))
    rep.append("H2 density fit rho = p / (%.6g + %.6g p), from %s.\n" % (Ah, Bh, os.path.basename(a.mix_table)))
    blocks_txt = "\n".join(b for n, b in pcs)
    mats_txt = "\n".join(b for n, b in keep)
    nppss = sum(1 for n, b in keep if n.startswith("ppss"))
    bnd_all = " ".join(aq["bnd"] for aq in aqs)

    aq_funcs, aq_fm, aq_bcs, aq_pp = [], [], [], []
    for aq in aqs:
        k = aq["k"]
        aq_funcs.append(f"""  [p_aq{k}_fn]   # Fetkovich aquifer {k}, re-based to the end of history: p0 = {aq['p0']:.0f} - {aq['cum_end']:.6g} / {aq['cap']}
    type = ParsedFunction
    symbol_names = 'p0 datum_elev rho_w g cum cap'
    symbol_values = '{aq['p0_new']:.2f} {aq['zd']} {aq['rho']} {aq['g']} aq{k}_cum_in_kg {aq['cap']}'
    expression = 'p0 - rho_w * g * (z - datum_elev) - cum / cap'
  []""")
        aq_fm.append(f"""  [aq{k}_influx]   # kg/s/m2 of water phase into the reservoir (+) or out (-); Eclipse's Fetkovich step average
    type = ADParsedFunctorMaterial
    property_name = aq{k}_influx
    functor_names = 'p_aq{k} pwater dt_size'
    functor_symbols = 'paq p dt'
    expression = '{aq['expr']}'
  []
  [aq{k}_water]   # component split of that flux: water ...
    type = ADParsedFunctorMaterial
    property_name = aq{k}_water
    functor_names = 'aq{k}_influx'
    functor_symbols = 'f'
    expression = '(1 - {EPS:g}) * f'
  []
  [aq{k}_ng]      # ... dissolved native gas (aquifer water carries native gas; outflow carries the cell's mix)
    type = ADParsedFunctorMaterial
    property_name = aq{k}_ng
    functor_names = 'aq{k}_influx zh2'
    functor_symbols = 'f h'
    expression = '{EPS:g} * f * if(f > 0, 1, 1 - min(max(h, 0), 1))'
  []
  [aq{k}_h2]      # ... dissolved H2 (only leaves the reservoir)
    type = ADParsedFunctorMaterial
    property_name = aq{k}_h2
    functor_names = 'aq{k}_influx zh2'
    functor_symbols = 'f h'
    expression = '{EPS:g} * f * if(f > 0, 0, min(max(h, 0), 1))'
  []""")
        for comp, var in (("water", "pwater"), ("ng", "sgas"), ("h2", "sgas"), ("h2", "zh2")):
            aq_bcs.append(f"""  [aquifer_{k}_{comp}_{var}]
    type = FVFunctorNeumannBC
    variable = {var}
    boundary = '{aq['bnd']}'
    functor = aq{k}_{comp}
  []""")
        aq_pp.append(f"""  [aq{k}_influx_kgps]   # + = water into the reservoir
    type = ADSideIntegralFunctorPostprocessor
    functor = aq{k}_influx
    boundary = '{aq['bnd']}'
    execute_on = 'initial timestep_end'
  []
  [aq{k}_cum_in_kg]
    type = TimeIntegratedPostprocessor
    value = aq{k}_influx_kgps
    time_integration_scheme = implicit-euler
    execute_on = 'initial timestep_end'
  []
  [aq{k}_ng_kgps]
    type = ADSideIntegralFunctorPostprocessor
    functor = aq{k}_ng
    boundary = '{aq['bnd']}'
    execute_on = 'initial timestep_end'
    outputs = none
  []
  [aq{k}_h2_kgps]
    type = ADSideIntegralFunctorPostprocessor
    functor = aq{k}_h2
    boundary = '{aq['bnd']}'
    execute_on = 'initial timestep_end'
    outputs = none
  []
  [aq{k}_cum_ng_in_kg]
    type = TimeIntegratedPostprocessor
    value = aq{k}_ng_kgps
    time_integration_scheme = implicit-euler
    execute_on = 'initial timestep_end'
  []
  [aq{k}_cum_h2_in_kg]
    type = TimeIntegratedPostprocessor
    value = aq{k}_h2_kgps
    time_integration_scheme = implicit-euler
    execute_on = 'initial timestep_end'
  []
  [aq{k}_pressure_at_datum]
    type = FunctionValuePostprocessor
    function = p_aq{k}_fn
    point = '0 0 {aq['zd']}'
    execute_on = 'initial timestep_end'
  []""")
    aq_aux = "\n".join(f"""  [p_aq{aq['k']}]
    type = MooseVariableFVReal
    two_term_boundary_expansion = false # aquifer pressure at the cell depth, as Eclipse
  []""" for aq in aqs)
    aq_auxk = "\n".join(f"""  [p_aq{aq['k']}]
    type = FunctionAux
    variable = p_aq{aq['k']}
    function = p_aq{aq['k']}_fn
    execute_on = 'initial timestep_begin'
  []""" for aq in aqs)
    sum_aq = lambda s: " + ".join(("aq%s_" % aq["k"]) + s for aq in aqs)
    aq_names = lambda s: " ".join(("aq%s_" % aq["k"]) + s for aq in aqs)
    na = len(aqs)
    sym = lambda p: " ".join("%s%d" % (p, i) for i in range(na))
    ssum = lambda p: " + ".join("%s%d" % (p, i) for i in range(na))
    ooc = (" -mat_mumps_icntl_22", " 1") if a.mumps_ooc else ("", "")

    nnc_mesh = nnc_fmat = nnc_k = ""
    if nnc_src:
        nnc_out = a.out + "_nnc.csv"
        if os.path.abspath(nnc_src) != os.path.abspath(nnc_out):
            shutil.copyfile(nnc_src, nnc_out)
        nnc_name = os.path.basename(nnc_out)
        rep.append("Non-neighbour (fault) connections of the history input carried over: %s "
                   "(FVPorousFlowNNCFlux; mobility lagged, pressure and composition implicit).\n" % nnc_name)
        nnc_mesh = "  allow_renumbering = false # the non-neighbour connections refer to the element ids of the file\n"
        nnc_fmat = f"""  # ---- composition on the non-neighbour (fault) connections: (water, native gas, H2) per phase ---
  [nnc_x_w1]
    type = ADParsedFunctorMaterial
    property_name = nnc_x_w1
    functor_names = 'zh2'
    functor_symbols = 'h'
    expression = '{EPS:g} * (1 - h)'
  []
  [nnc_x_w2]
    type = ADParsedFunctorMaterial
    property_name = nnc_x_w2
    functor_names = 'zh2'
    functor_symbols = 'h'
    expression = '{EPS:g} * h'
  []
  [nnc_x_g1]
    type = ADParsedFunctorMaterial
    property_name = nnc_x_g1
    functor_names = 'zh2'
    functor_symbols = 'h'
    expression = '1 - h'
  []
"""
        common = (f"    nnc_file = {nnc_name}\n    pressure = pwater\n    pressure_offsets = '0 pc_lag'\n"
                  "    relperm = 'krw krg'\n    density = 'water_density rho_g_lag'\n    viscosity = 'muw mug'\n")
        for nm, var, xf, cm in (("nnc_water", "pwater", "'%g 0'" % (1 - EPS), "water"),
                                ("nnc_ng", "sgas", "'nnc_x_w1 nnc_x_g1'", "native gas, total-gas equation"),
                                ("nnc_h2_tot", "sgas", "'nnc_x_w2 zh2'", "H2, total-gas equation"),
                                ("nnc_h2", "zh2", "'nnc_x_w2 zh2'", "H2 equation")):
            nnc_k += (f"  [{nm}]   # {cm}: flow through the non-neighbour (fault) connections\n"
                      f"    type = FVPorousFlowNNCFlux\n    variable = {var}\n{common}    mass_fraction = {xf}\n  []\n")
    hyst_auxv = hyst_auxk = ""
    if hyst:
        hyst_auxv = """  [sg_max]   # largest gas saturation each cell has had (rel-perm hysteresis), from the history run
    type = MooseVariableFVReal
    initial_from_file_var = sg_max
  []
"""
        hyst_auxk = """  [sg_max]   # updated once per time step (lagged); the gas rel-perm follows Killough scanning curves below it
    type = ParsedAux
    variable = sg_max
    coupled_variables = 'sg_max sgas'
    expression = 'max(sg_max, sgas)'
    execute_on = timestep_end
  []
"""
    mix_zh2 = os.path.splitext(os.path.basename(a.mix_table))[0] + "_zh2.csv"
    s = f"""# =============================================================================================
# Cyclic hydrogen storage, finite-volume PorousFlow model (two-point flux, as the history match)
#   {a.cycles} cycles of: inject {a.mass_t:g} t of pure H2 over {a.inj_days:g} d, then withdraw {a.mass_t:g} t of
#   the gas mixture over {a.wd_days:g} d.  Well {well}.  Start = end of the history run
#   {os.path.basename(a.hist_exodus)} (t = {times[step]:.0f} s): pressure and saturation of every cell.
#   Written by gen_h2_fv.py from {os.path.basename(a.hist_input)}; grid, rock, rel-perm, capillary pressure,
#   water PVT and aquifer are copied from that input.
#
#   3 components, 2 phases:
#     component 0 = water (equation on pwater), 2 = hydrogen (on zh2); total gas 1 + 2 on sgas
#     phase 0 = water (+ {EPS:g} dissolved gas, see [massfrac]); phase 1 = gas (native gas + H2)
#   primary variables (one value per cell): pwater [Pa], sgas [-], zh2 = H2 mass fraction of the gas
#   (xng = 1 - zh2, the native-gas fraction, is kept as an output variable)
#   Mixture density and viscosity rho(p, zh2), mu(p, zh2) at the reservoir temperature are tabulated in
#   {mix_zh2}; the table interpolates in (p, T), so zh2 is passed in the temperature slot
#   (the model is isothermal; nothing else uses temperature: water density has no thermal expansion).
#   Mixing: advection of the mixture plus the numerical dispersion of the upwind two-point scheme
#   (about half a cell, ~20 m laterally, ~0.5 m vertically on the full grid), larger than molecular
#   diffusion (~3 m in 6 months); no explicit dispersion kernels.
#   Well: target rates, limited by the bottom-hole pressure (BHP): during injection the BHP may not exceed
#   {p_max / PSI_:.1f} psia (the initial reservoir pressure at the BHP reference depth), during withdrawal it
#   may not fall below {p_min / PSI_:.1f} psia. Withdrawn gas comes from each completion cell and is split
#   between native gas and H2 by that cell's composition (implicit), so purity and recovery are results.
#   Requires an app with the add-ons FVPorousFlowAdvectiveFluxTPFA, ADPorousFlowPorosityRock,
#   ADPorousFlowMassFractionLinear and FVMaxIncrementDamper (and FVPorousFlowNNCFlux with fault NNCs).
# =============================================================================================

[Mesh]
  [file]
    type = FileMeshGenerator
    file = {os.path.basename(out_e)}
    use_for_exodus_restart = true
  []
{nnc_mesh}[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 -9.81'
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pwater sgas zh2'
    number_fluid_phases = 2
    number_fluid_components = 3
  []
{blocks_txt}
[]

[Variables]
  [pwater]
    type = MooseVariableFVReal
    two_term_boundary_expansion = false # aquifer flux uses the cell pressure, as Eclipse
    initial_from_file_var = p_init # end-of-history water pressure
  []
  [sgas]
    type = MooseVariableFVReal
    initial_from_file_var = sg_init # end-of-history gas saturation
  []
  [zh2]   # H2 mass fraction of the gas
    type = MooseVariableFVReal
    two_term_boundary_expansion = false
    initial_condition = 0 # only native gas at the start
  []
[]

[Functions]
{chr(10).join(aq_funcs)}
  [inj_fn]    # H2 injection rate [kg/s]; evaluated 1 s before the step end so a step ending on a switch belongs to its period
    type = ParsedFunction
    symbol_names = 'P Ti q tend'
    symbol_values = '{P:.1f} {Ti:.1f} {q_inj:.9g} {t_end:.1f}'
    expression = 'if((t - 1) - floor((t - 1) / P) * P < Ti & t <= tend, q, 0)'
  []
  [wd_fn]     # gas-mixture withdrawal rate [kg/s, positive]
    type = ParsedFunction
    symbol_names = 'P Ti q'
    symbol_values = '{P:.1f} {Ti:.1f} {q_wd:.9g}'
    expression = 'if(t > 0 & (t - 1) - floor((t - 1) / P) * P >= Ti, q, 0)'
  []
  [schedule]  # switch times: the time stepper lands on them
    type = PiecewiseLinear
    x = '0 {sw_s}'
    y = '{sched_y}'
  []
[]

[AuxVariables]
{hyst_auxv}  [bounds_dummy]   # needed by the [Bounds] objects
    type = MooseVariableFVReal
  []
  [one]
    type = MooseVariableFVReal
    initial_condition = 1
  []
  [zero]
    type = MooseVariableFVReal
    initial_condition = 0
  []
  [poro]
    type = MooseVariableFVReal
    initial_from_file_var = poro
  []
  [perm_xx]
    type = MooseVariableFVReal
    initial_from_file_var = perm_xx
  []
  [perm_yy]
    type = MooseVariableFVReal
    initial_from_file_var = perm_yy
  []
  [perm_zz]
    type = MooseVariableFVReal
    initial_from_file_var = perm_zz
  []
  [perm_xy]
    type = MooseVariableFVReal
    initial_from_file_var = perm_xy
  []
  [ww_{well}]   # rate fraction / cell volume in the completion cells [1/m3]
    type = MooseVariableFVReal
    initial_from_file_var = ww_{well}
  []
{aq_aux}
  [pgas]
    type = MooseVariableFVReal
  []
  [xh2]       # H2 mass fraction of the gas (0 where there is no gas)
    type = MooseVariableFVReal
  []
  [xng]       # native-gas mass fraction of the gas, 1 - zh2 (output only)
    type = MooseVariableFVReal
  []
  [yh2]       # H2 mole fraction of the gas
    type = MooseVariableFVReal
  []
  [gas_density]
    type = MooseVariableFVReal
  []
  [krg]       # well mobility pieces, from the last converged step
    type = MooseVariableFVReal
  []
  [krw]
    type = MooseVariableFVReal
  []
  [mug]
    type = MooseVariableFVReal
    initial_condition = 1.5e-5
  []
  [muw]
    type = MooseVariableFVReal
    initial_condition = 2.8e-4
  []
  [pc_lag]    # capillary pressure pgas - pwater
    type = MooseVariableFVReal
  []
  [rho_g_lag] # gas density
    type = MooseVariableFVReal
  []
  [bhp_inj_cell]
    type = MooseVariableFVReal
  []
  [bhp_wd_cell]
    type = MooseVariableFVReal
  []
  [water_density]
    type = MooseVariableFVReal
  []
  [phi]       # porosity (with rock compressibility)
    type = MooseVariableFVReal
  []
  [h2_md]     # kg of H2 per m3 of rock
    type = MooseVariableFVReal
  []
  [ng_md]     # kg of native gas per m3 of rock
    type = MooseVariableFVReal
  []
  [w_md]      # kg of water per m3 of rock
    type = MooseVariableFVReal
  []
  [hcpv]      # gas pore volume per m3 of rock
    type = MooseVariableFVReal
  []
  [p_hcpv]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
{aq_auxk}{hyst_auxk}  [pgas]
    type = ADPorousFlowPropertyAux
    variable = pgas
    property = pressure
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [gas_density]
    type = ADPorousFlowPropertyAux
    variable = gas_density
    property = density
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [krg]
    type = ADPorousFlowPropertyAux
    variable = krg
    property = relperm
    phase = 1
    execute_on = 'initial timestep_begin' # lagged: the same values during the step and at its end
  []
  [krw]
    type = ADPorousFlowPropertyAux
    variable = krw
    property = relperm
    phase = 0
    execute_on = 'initial timestep_begin' # lagged: the same values during the step and at its end
  []
  [mug]
    type = ADPorousFlowPropertyAux
    variable = mug
    property = viscosity
    phase = 1
    execute_on = 'initial timestep_begin' # lagged: the same values during the step and at its end
  []
  [muw]
    type = ADPorousFlowPropertyAux
    variable = muw
    property = viscosity
    phase = 0
    execute_on = 'initial timestep_begin' # lagged: the same values during the step and at its end
  []
  [rho_g_lag]
    type = ADPorousFlowPropertyAux
    variable = rho_g_lag
    property = density
    phase = 1
    execute_on = 'initial timestep_begin'
  []
  [pc_lag]
    type = ParsedAux
    variable = pc_lag
    coupled_variables = 'pgas pwater'
    expression = 'pgas - pwater'
    execute_on = 'initial timestep_begin' # lagged: the same values during the step and at its end
  []
  [bhp_inj_cell]
    type = ADFunctorElementalAux
    variable = bhp_inj_cell
    functor = bhp_inj
    execute_on = 'timestep_end'
  []
  [bhp_wd_cell]
    type = ADFunctorElementalAux
    variable = bhp_wd_cell
    functor = bhp_wd
    execute_on = 'timestep_end'
  []
  [water_density]
    type = ADPorousFlowPropertyAux
    variable = water_density
    property = density
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [phi]
    type = ADPorousFlowPropertyAux
    variable = phi
    property = porosity
    execute_on = 'initial timestep_end'
  []
  [xng]
    type = ParsedAux
    variable = xng
    coupled_variables = 'zh2'
    expression = '1 - zh2'
    execute_on = 'initial timestep_end'
  []
  [xh2]
    type = ParsedAux
    variable = xh2
    coupled_variables = 'zh2 sgas'
    expression = 'if(sgas > 1e-6, zh2, 0)'
    execute_on = 'initial timestep_end'
  []
  [yh2]
    type = ParsedAux
    variable = yh2
    coupled_variables = 'xh2'
    expression = 'xc := min(max(xh2, 0), 1); (xc / {M_H2}) / (xc / {M_H2} + (1 - xc) / {MNG})'
    execute_on = 'initial timestep_end'
  []
  [h2_md]
    type = ParsedAux
    variable = h2_md
    coupled_variables = 'phi sgas zh2 gas_density water_density'
    expression = 'phi * (sgas * gas_density * zh2 + (1 - sgas) * water_density * {EPS:g} * zh2)'
    execute_on = 'initial timestep_end'
  []
  [ng_md]
    type = ParsedAux
    variable = ng_md
    coupled_variables = 'phi sgas zh2 gas_density water_density'
    expression = 'phi * (sgas * gas_density * (1 - zh2) + (1 - sgas) * water_density * {EPS:g} * (1 - zh2))'
    execute_on = 'initial timestep_end'
  []
  [w_md]
    type = ParsedAux
    variable = w_md
    coupled_variables = 'phi sgas water_density'
    expression = 'phi * (1 - sgas) * water_density * (1 - {EPS:g})'
    execute_on = 'initial timestep_end'
  []
  [hcpv]
    type = ParsedAux
    variable = hcpv
    coupled_variables = 'phi sgas'
    expression = 'phi * max(sgas, 0)'
    execute_on = 'initial timestep_end'
  []
  [p_hcpv]
    type = ParsedAux
    variable = p_hcpv
    coupled_variables = 'phi sgas pgas'
    expression = 'phi * max(sgas, 0) * pgas'
    execute_on = 'initial timestep_end'
  []
[]

[FVKernels]
  [mass_water_dot]
    type = FVPorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pwater
  []
  [flux_water]
    type = FVPorousFlowAdvectiveFluxTPFA
    fluid_component = 0
    variable = pwater
  []
  # Equations: water on pwater; TOTAL gas (native gas + H2) on sgas; H2 on zh2. Each conservation law
  # is kept exactly (the native-gas one is total minus H2). The total-gas equation always depends on sgas,
  # and the H2 equation increases with zh2, as the bounded solver needs: when zh2 sits on a bound (0: no
  # H2, 1: no native gas) the solver drops the H2 equation and the total-gas equation, which is then the
  # conservation law of the one gas present, keeps sgas determined.
  [mass_ng_dot]   # native gas, in the total-gas equation (on sgas)
    type = FVPorousFlowMassTimeDerivative
    fluid_component = 1
    variable = sgas
  []
  [flux_ng]
    type = FVPorousFlowAdvectiveFluxTPFA
    fluid_component = 1
    variable = sgas
  []
  [mass_h2_dot_tot]   # H2, in the total-gas equation (on sgas)
    type = FVPorousFlowMassTimeDerivative
    fluid_component = 2
    variable = sgas
  []
  [flux_h2_tot]
    type = FVPorousFlowAdvectiveFluxTPFA
    fluid_component = 2
    variable = sgas
  []
  [mass_h2_dot]   # H2 equation (on zh2)
    type = FVPorousFlowMassTimeDerivative
    fluid_component = 2
    variable = zh2
  []
  [flux_h2]
    type = FVPorousFlowAdvectiveFluxTPFA
    fluid_component = 2
    variable = zh2
  []
{nnc_k}  [h2_injection]
    type = FVCoupledForce
    variable = zh2
    v = src_h2_inj
  []
  [h2_injection_tot]
    type = FVCoupledForce
    variable = sgas
    v = src_h2_inj
  []
  [h2_withdrawal]
    type = FVCoupledForce
    variable = zh2
    v = src_h2_wd
  []
  [h2_withdrawal_tot]
    type = FVCoupledForce
    variable = sgas
    v = src_h2_wd
  []
  [ng_withdrawal]
    type = FVCoupledForce
    variable = sgas
    v = src_ng_wd
  []
[]

[FVBCs]
{chr(10).join(aq_bcs)}
[]

[FluidProperties]
{water}
  [gas]   # native gas + H2 mixture; 'temperature' column = H2 mass fraction zh2
    type = TabulatedBicubicFluidProperties
    fluid_property_file = {mix_zh2}
    interpolated_properties = 'density viscosity'
    out_of_bounds_behavior = set_to_closest_bound
  []
[]

[Materials]
  [temperature]   # composition carried in the temperature slot (see header)
    type = ADPorousFlowTemperature
    temperature = zh2
  []
  [massfrac]
    # gas phase:   (water, native gas, H2) = (0, 1 - zh2, zh2)
    # water phase: (water, native gas, H2) = (1 - EPS, EPS (1 - zh2), EPS zh2), EPS = {EPS:g}
    # The small dissolved fraction gives zh2 a conservation equation where there is no free gas, so
    # Newton stays well posed at the gas front. (Real solubilities: CH4 ~1.5e-3, H2 ~2e-6 by mass.)
    type = ADPorousFlowMassFractionLinear
    mass_fraction_vars = 'one zh2 zero zh2'
    multipliers = '1 {-EPS:g} 1 -1'
    offsets = '{-EPS:g} {EPS:g} 0 1'
  []
  [water]
    type = ADPorousFlowSingleComponentFluid
    fp = water
    phase = 0
    compute_enthalpy = false
    compute_internal_energy = false
  []
  [gas]
    type = ADPorousFlowSingleComponentFluid
    fp = gas
    phase = 1
    compute_enthalpy = false
    compute_internal_energy = false
  []
{mats_txt}
[]

[FunctorMaterials]
{nnc_fmat}  # ---- well: rate target, limited by the bottom-hole pressure ---------------------------------------
  # Each completion cell c takes its share w_c V_c of the target rate (w = ww_{well} = share / cell volume,
  # shares from COMPDAT KH, proportional to the connection factors CF), unless that would need a well
  # pressure beyond the limit; then the cell gets the rate the limiting pressure gives:
  #   q_c = CF_c * rho * mobility * |p_limit(z_c) - p_gas,c|,   CF_c = w_c V_c * sum(CF) = w_c V_c * {C_well:.6g} m3
  # p_limit(z) is the BHP limit at the reference depth z_ref = {-z_ref_m:.2f} m moved along the initial gas
  # gradient ({rho_ref:.1f} kg/m3): p_max(z) = {p_max:.0f} + {rho_ref:.1f} g (z_ref - z),  p_min(z) likewise with {p_min:.0f}.
  # Mobility (kr/mu) is taken from the last converged step; pressures and composition are implicit.
  [pg_cell]
    type = ADParsedFunctorMaterial
    property_name = pg_cell
    functor_names = 'pwater pc_lag'
    functor_symbols = 'p pc'
    expression = 'p + pc'
  []
  [rho_h2]       # density of the injected H2 at the cell's gas pressure
    type = ADParsedFunctorMaterial
    property_name = rho_h2
    functor_names = 'pg_cell'
    functor_symbols = 'p'
    expression = 'p / ({Ah:.8g} + {Bh:.8g} * p)'
  []
  [inj_cap]      # kg/s/m3 the cell could take with the well at the injection limit
    type = ADParsedFunctorMaterial
    property_name = inj_cap
    functor_names = 'ww_{well} rho_h2 krg mug krw muw pg_cell'
    functor_symbols = 'w rho krg mug krw muw p'
    expression = 'w * {C_well:.8g} * rho * (krg / mug + krw / muw) * max({p_max:.2f} + {rho_ref:.4f} * 9.81 * ({-z_ref_m:.4f} - z) - p, 0)'
  []
  [wd_cap]       # kg/s/m3 of gas the cell could deliver with the well at the withdrawal limit (deliverability)
    type = ADParsedFunctorMaterial
    property_name = wd_cap
    functor_names = 'ww_{well} rho_g_lag krg mug pg_cell'
    functor_symbols = 'w rho krg mug p'
    expression = 'w * {C_well:.8g} * rho * krg / mug * max(p - {p_min:.2f} - {rho_ref:.4f} * 9.81 * ({-z_ref_m:.4f} - z), 0)'
  []
  [src_h2_inj]   # kg/s/m3 of H2 injected
    type = ADParsedFunctorMaterial
    property_name = src_h2_inj
    functor_names = 'ww_{well} inj_fn inj_cap'
    functor_symbols = 'w q cap'
    expression = 'min(w * q, cap)'
  []
  [wd_cell]      # kg/s/m3 of gas withdrawn (positive)
    type = ADParsedFunctorMaterial
    property_name = wd_cell
    functor_names = 'ww_{well} wd_fn wd_cap'
    functor_symbols = 'w q cap'
    expression = 'min(w * q, cap)'
  []
  [src_h2_wd]    # kg/s/m3 of H2 withdrawn (negative): the cell's withdrawal x its H2 fraction
    type = ADParsedFunctorMaterial
    property_name = src_h2_wd
    functor_names = 'wd_cell zh2'
    functor_symbols = 'qc h'
    expression = '-qc * min(max(h, 0), 1)'
  []
  [src_ng_wd]    # kg/s/m3 of native gas withdrawn (negative)
    type = ADParsedFunctorMaterial
    property_name = src_ng_wd
    functor_names = 'wd_cell zh2'
    functor_symbols = 'qc h'
    expression = '-qc * (1 - min(max(h, 0), 1))'
  []
  [deliv_sm3pd]  # deliverability density: wd_cap in standard m3 per day (mix of the cell's gas)
    type = ADParsedFunctorMaterial
    property_name = deliv_sm3pd
    functor_names = 'wd_cap zh2'
    functor_symbols = 'cap h'
    expression = 'cap * 86400 * ((1 - min(max(h, 0), 1)) / {rho_sc_ng:.6g} + min(max(h, 0), 1) / {a.rho_sc_h2:.6g})'
  []
  [bhp_inj]      # BHP at the reference depth implied by each completion cell while injecting (-1e30 elsewhere)
    type = ADParsedFunctorMaterial
    property_name = bhp_inj
    functor_names = 'ww_{well} src_h2_inj rho_h2 krg mug krw muw pg_cell'
    functor_symbols = 'w q rho krg mug krw muw p'
    expression = 'if(w > 0, p + q / (w * {C_well:.8g} * rho * (krg / mug + krw / muw) + 1e-30) + {rho_ref:.4f} * 9.81 * (z - {-z_ref_m:.4f}), -1e30)'
  []
  [bhp_wd]       # same while withdrawing (+1e30 elsewhere)
    type = ADParsedFunctorMaterial
    property_name = bhp_wd
    functor_names = 'ww_{well} wd_cell rho_g_lag krg mug pg_cell'
    functor_symbols = 'w q rho krg mug p'
    expression = 'if(w > 0, p - q / (w * {C_well:.8g} * rho * krg / mug + 1e-30) + {rho_ref:.4f} * 9.81 * (z - {-z_ref_m:.4f}), 1e30)'
  []
{chr(10).join(aq_fm)}
[]

[Postprocessors]
  [dt_size]
    type = TimestepSize
    execute_on = 'initial timestep_begin'
  []
  # ---- well ------------------------------------------------------------------------------------
  [h2_inj_target_kgps]
    type = FunctionValuePostprocessor
    function = inj_fn
    execute_on = 'initial timestep_end'
  []
  [gas_wd_target_kgps]
    type = FunctionValuePostprocessor
    function = wd_fn
    execute_on = 'initial timestep_end'
  []
  [h2_inj_rate_kgps]   # actual (pressure-limited) injection rate
    type = ADElementIntegralFunctorPostprocessor
    functor = src_h2_inj
    execute_on = 'initial timestep_end'
  []
  [gas_wd_rate_kgps]   # actual (pressure-limited) withdrawal rate
    type = ADElementIntegralFunctorPostprocessor
    functor = wd_cell
    execute_on = 'initial timestep_end'
  []
  [deliverability_sm3pd]   # gas rate the well could deliver now at the minimum BHP
    type = ADElementIntegralFunctorPostprocessor
    functor = deliv_sm3pd
    execute_on = 'initial timestep_end'
  []
  [injectivity_sm3pd]   # H2 rate the well could take now with the BHP at the injection limit (standard m3/d)
    type = ADElementIntegralFunctorPostprocessor
    functor = inj_cap
    prefactor = {86400.0 / a.rho_sc_h2:.6g}
    execute_on = 'initial timestep_end'
  []
  [withdrawal_sm3pd]   # actual withdrawal rate in standard m3/day
    type = ParsedPostprocessor
    pp_names = 'h2_wd_rate_kgps ng_wd_rate_kgps'
    pp_symbols = 'h n'
    expression = '86400 * (h / {a.rho_sc_h2:.6g} + n / {rho_sc_ng:.6g})'
    execute_on = 'initial timestep_end'
  []
  [bhp_inj_max]
    type = ElementExtremeValue
    variable = bhp_inj_cell
    value_type = max
    execute_on = 'timestep_end'
    outputs = none
  []
  [bhp_wd_min]
    type = ElementExtremeValue
    variable = bhp_wd_cell
    value_type = min
    execute_on = 'timestep_end'
    outputs = none
  []
  [bhp_Pa]   # bottom-hole pressure at the reference depth (injection: highest, withdrawal: lowest over the completions)
    type = ParsedPostprocessor
    pp_names = 'bhp_inj_max bhp_wd_min h2_inj_target_kgps'
    pp_symbols = 'bi bw qi'
    expression = 'if(qi > 0, bi, bw)'
    execute_on = 'timestep_end'
  []
  [bhp_max_limit_Pa]
    type = ConstantPostprocessor
    value = {p_max:.2f}
    execute_on = 'initial timestep_end'
  []
  [bhp_min_limit_Pa]
    type = ConstantPostprocessor
    value = {p_min:.2f}
    execute_on = 'initial timestep_end'
  []
  [h2_wd_rate_kgps]   # H2 withdrawal rate (positive)
    type = ADElementIntegralFunctorPostprocessor
    functor = src_h2_wd
    prefactor = -1
    execute_on = 'initial timestep_end'
  []
  [ng_wd_rate_kgps]   # native-gas withdrawal rate (positive)
    type = ADElementIntegralFunctorPostprocessor
    functor = src_ng_wd
    prefactor = -1
    execute_on = 'initial timestep_end'
  []
  [cum_h2_injected_kg]
    type = TimeIntegratedPostprocessor
    value = h2_inj_rate_kgps
    time_integration_scheme = implicit-euler
    execute_on = 'initial timestep_end'
  []
  [cum_h2_withdrawn_kg]
    type = TimeIntegratedPostprocessor
    value = h2_wd_rate_kgps
    time_integration_scheme = implicit-euler
    execute_on = 'initial timestep_end'
  []
  [cum_ng_withdrawn_kg]
    type = TimeIntegratedPostprocessor
    value = ng_wd_rate_kgps
    time_integration_scheme = implicit-euler
    execute_on = 'initial timestep_end'
  []
  [h2_recovery_cum]   # cumulative H2 withdrawn / cumulative H2 injected
    type = ParsedPostprocessor
    pp_names = 'cum_h2_withdrawn_kg cum_h2_injected_kg'
    pp_symbols = 'w i'
    expression = 'if(i > 0, w / i, 0)'
    execute_on = 'initial timestep_end'
  []
  [h2_molefrac_withdrawn]   # H2 mole fraction of the gas being withdrawn
    type = ParsedPostprocessor
    pp_names = 'h2_wd_rate_kgps ng_wd_rate_kgps'
    pp_symbols = 'h n'
    expression = 'if(h + n > 0, (h / {M_H2}) / (h / {M_H2} + n / {MNG}), 0)'
    execute_on = 'initial timestep_end'
  []
  [xng_well]   # native-gas mass fraction in the top completion cell
    type = PointValue
    variable = xng
    point = '{pt}'
    execute_on = 'initial timestep_end'
  []
  [p_well_Pa]   # gas pressure in the top completion cell
    type = PointValue
    variable = pgas
    point = '{pt}'
    execute_on = 'initial timestep_end'
  []
  [p_well_max_Pa]
    type = TimeExtremeValue
    postprocessor = p_well_Pa
    value_type = max
    execute_on = 'initial timestep_end'
  []
  # ---- reservoir -------------------------------------------------------------------------------
  [hcpv_m3]
    type = ElementIntegralVariablePostprocessor
    variable = hcpv
    execute_on = 'initial timestep_end'
  []
  [p_hcpv_int]
    type = ElementIntegralVariablePostprocessor
    variable = p_hcpv
    execute_on = 'initial timestep_end'
    outputs = none
  []
  [p_hcpv_Pa]   # gas-pore-volume-weighted gas pressure (as Eclipse FPR)
    type = ParsedPostprocessor
    pp_names = 'p_hcpv_int hcpv_m3'
    pp_symbols = 'a b'
    expression = 'a / b'
    execute_on = 'initial timestep_end'
  []
{chr(10).join(aq_pp)}
  # ---- inventories and mass balance (all from the 3D fields; must equal the well/aquifer totals) --
  [h2_in_place_kg]
    type = ElementIntegralVariablePostprocessor
    variable = h2_md
    execute_on = 'initial timestep_end'
  []
  [ng_in_place_kg]
    type = ElementIntegralVariablePostprocessor
    variable = ng_md
    execute_on = 'initial timestep_end'
  []
  [water_in_place_kg]
    type = ElementIntegralVariablePostprocessor
    variable = w_md
    execute_on = 'initial timestep_end'
  []
  [ng_in_place_initial_kg]
    type = ElementIntegralVariablePostprocessor
    variable = ng_md
    execute_on = 'initial'
    outputs = none
  []
  [water_in_place_initial_kg]
    type = ElementIntegralVariablePostprocessor
    variable = w_md
    execute_on = 'initial'
    outputs = none
  []
  [h2_balance_kg]   # H2 in place - (injected - withdrawn + aquifer): should stay ~0
    type = ParsedPostprocessor
    pp_names = 'h2_in_place_kg cum_h2_injected_kg cum_h2_withdrawn_kg {aq_names('cum_h2_in_kg')}'
    pp_symbols = 'm i w {sym('a')}'
    expression = 'm - (i - w + {ssum('a')})'
    execute_on = 'timestep_end'
  []
  [ng_balance_kg]
    type = ParsedPostprocessor
    pp_names = 'ng_in_place_kg ng_in_place_initial_kg cum_ng_withdrawn_kg {aq_names('cum_ng_in_kg')}'
    pp_symbols = 'm m0 w {sym('a')}'
    expression = 'm - (m0 - w + {ssum('a')})'
    execute_on = 'timestep_end'
  []
  [water_balance_kg]
    type = ParsedPostprocessor
    pp_names = 'water_in_place_kg water_in_place_initial_kg {aq_names('cum_in_kg')}'
    pp_symbols = 'm m0 {sym('a')}'
    expression = 'm - (m0 + (1 - {EPS:g}) * ({ssum('a')}))'
    execute_on = 'timestep_end'
  []
  [min_xng]
    type = ElementExtremeValue
    variable = xng
    value_type = min
    execute_on = 'initial timestep_end'
  []
  [max_xng]
    type = ElementExtremeValue
    variable = xng
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [nl_its]
    type = NumNonlinearIterations
  []
[]

[Bounds]   # 0 <= zh2 <= 1 in every cell (no negative H2, no negative native gas), enforced by PETSc's
           # variational-inequality Newton solver (vinewtonrsls). Cells with almost no gas determine zh2 only
           # through the tiny dissolved fraction; without the bounds Newton can push it outside [0, 1] and the
           # upwind fluxes then spread negative component masses. On a bound the H2 equation is dropped and
           # the total-gas equation (on sgas) carries the conservation of the gas that is present.
  [zh2_lo]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = zh2
    bound_type = lower
    bound_value = 0
  []
  [zh2_hi]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = zh2
    bound_type = upper
    bound_value = 1
  []
[]

[Dampers]   # Newton damping (add-on FVMaxIncrementDamper; MOOSE's MaxIncrement is for finite-element variables)
  [limit_p]
    type = FVMaxIncrementDamper
    variable = pwater
    max_increment = 1e6 # Pa per Newton iteration
  []
  [limit_s]
    type = FVMaxIncrementDamper
    variable = sgas
    max_increment = 0.2
  []
  [limit_x]   # composition: at most 0.2 per iteration where there is gas ([0, 1] is enforced in every cell by [Bounds])
    type = FVMaxIncrementDamper
    variable = zh2
    max_increment = 0.2
    ignore_variable = sgas # composition is poorly determined in cells with almost no gas: not damped there
    ignore_below = 1e-3
  []
[]

[Preconditioning]
  [mumps]   # direct solver; ~3 unknowns per cell
    type = SMP
    full = true
    petsc_options_iname = '-snes_type    -pc_type -pc_factor_mat_solver_type -mat_mumps_icntl_14{ooc[0]}'
    petsc_options_value = ' vinewtonrsls  lu       mumps                      {a.mumps_relax}{ooc[1]}'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = none
  automatic_scaling = false
  end_time = {t_end:.1f} # s ({a.cycles} cycles x {(a.inj_days + a.wd_days):g} d = {t_end / DAY:g} d)
  dtmax = {a.dtmax_days * DAY:.0f} # {a.dtmax_days:g} days
  dtmin = 1 # s: stop (checkpoint kept) instead of cutting the step forever if Newton cannot converge
  nl_rel_tol = {a.nl_rel_tol:g}
  nl_abs_tol = {a.nl_abs_tol:g} # kg/s summed over the cells
  nl_max_its = 12
  l_max_its = 100
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 3600
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 6
    iteration_window = 2
    timestep_limiting_function = schedule
    force_step_every_function_point = true
  []
[]

[Outputs]
  checkpoint = true
  print_linear_residuals = false
  [csv]
    type = CSV
  []
  [exodus]   # 3D fields at the start, at the end of every injection and withdrawal period{"" if a.exodus_days <= 0 else " and every %g days" % a.exodus_days}
    type = Exodus
    sync_times = '{ex_s}'
    sync_only = true
    execute_on = 'initial timestep_end'
    hide = 'one zero bounds_dummy poro perm_xx perm_yy perm_zz perm_xy ww_{well} h2_md ng_md w_md hcpv p_hcpv krg krw mug muw pc_lag rho_g_lag bhp_inj_cell bhp_wd_cell {" ".join("p_aq%s" % aq["k"] for aq in aqs)}'
  []
[]
"""
    open(a.out + ".i", "w", newline="\n").write(s)
    write_zh2_table(a.mix_table, os.path.join(os.path.dirname(os.path.abspath(a.out)), mix_zh2))
    rep.append("Schedule: %d cycles, injection %.6g kg/s for %g d, withdrawal %.6g kg/s for %g d; end %.0f s.\n"
               % (a.cycles, q_inj, a.inj_days, q_wd, a.wd_days, t_end))
    rep.append("Copied from the history input: %d capillary-pressure objects, %d ppss materials, %d other rock materials; "
               "well %s (point %s); aquifer boundaries %s.\n" % (len(pcs), nppss, len(keep) - nppss, well, pt, bnd_all))
    open(a.out + "_report.md", "w", newline="\n").write("\n".join(rep))
    print("".join(rep))
    print("wrote %s.i, %s.e, %s_report.md" % (a.out, a.out, a.out))


if __name__ == "__main__":
    main()

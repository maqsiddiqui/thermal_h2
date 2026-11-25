#!/usr/bin/env python3
"""compare_with_eclipse.py - compare a MOOSE run made with em2ex2moose.py against the Eclipse results.

From the Miniforge prompt:
  python compare_with_eclipse.py --eclipse-h5 CROFT_files/HM_MODEL.h5 --run "MOOSE (em2ex mesh)=croft_em2ex"

  --run LABEL=PREFIX   a MOOSE run; PREFIX is the -o prefix given to em2ex2moose.py. The script reads
                       PREFIX_csv.csv (time series), PREFIX_meta.json (units, densities), and, for the
                       field pressure, PREFIX_exodus.e (MOOSE results) with PREFIX.e (mesh).
                       Repeat --run to put several runs on the same plots.
  -o FILE.png          where to save the figure (default moose_vs_eclipse.png)
  --no-field-pressure  skip the field-pressure curve (reading PREFIX_exodus.e takes a minute or two)
  --days 0 180 365     days for the printed table

In a Jupyter notebook, use the functions directly (see compare_with_eclipse.ipynb):
  import compare_with_eclipse as cwe
  ecl  = cwe.read_eclipse("CROFT_files/HM_MODEL.h5", units="FIELD")
  runs = {"MOOSE": cwe.read_moose("croft_em2ex")}
  fig, ax = cwe.plot_comparison(ecl, runs)
  print(cwe.summary_table(ecl, runs))

Everything is converted to field units for the plots: Bscf, MMscf/day, psia, and million stb for water
(water in place at stock-tank conditions, as Eclipse FWIP).
"""
from __future__ import annotations

import argparse
import csv
import json
import os

import numpy as np

PSI = 6894.757293168          # Pa per psi
BAR = 1.0e5                   # Pa per bar
SCF = 0.028316846592          # m3 per scf
STB = 0.158987294928          # m3 per stb / rb

# Eclipse summary units -> SI, per unit system
ECL_UNITS = {
    "FIELD": dict(gas_vol=1000 * SCF, liq_vol=STB, pressure=PSI),      # Mscf, stb, psia
    "METRIC": dict(gas_vol=1.0, liq_vol=1.0, pressure=BAR),            # sm3, sm3, barsa
}


# ----------------------------------------------------------------------------------------------
# reading
# ----------------------------------------------------------------------------------------------
def read_eclipse(h5_path, units="FIELD", well=None):
    """Eclipse summary from the .h5 file Petrel writes, converted to SI.

    Returns a dict of numpy arrays: t_days, gip_sm3, gas_cum_sm3, gas_rate_sm3d, p_field_pa,
    p_well_block_pa, bhp_pa, wip_sm3 (water in place, stock-tank m3). Missing vectors are left out."""
    import h5py
    u = ECL_UNITS[units.upper()]
    out = {}
    with h5py.File(h5_path, "r") as f:
        out["t_days"] = f["general/time"][:].astype(float)
        sv = f["summary_vectors"]

        def first(key, idx=0):
            if key not in sv:
                return None
            ids = sorted(sv[key].keys(), key=lambda s: int(s) if s.isdigit() else s)
            return sv[key][ids[idx]]["values"][:].astype(float) if ids else None

        pairs = [("gip_sm3", "FGIP", "gas_vol"), ("gas_cum_sm3", "FGPT", "gas_vol"),
                 ("gas_rate_sm3d", "FGPR", "gas_vol"), ("p_field_pa", "FPR", "pressure"),
                 ("p_well_block_pa", "WBP", "pressure"), ("bhp_pa", "WBHP", "pressure"),
                 ("wip_sm3", "FWIP", "liq_vol"), ("water_cum_sm3", "FWPT", "liq_vol"),
                 ("water_rate_sm3d", "WWPR", "liq_vol"), ("water_prod_sm3", "WWPT", "liq_vol")]
        for name, key, kind in pairs:
            v = first(key)
            if v is not None:
                out[name] = v * u[kind]
    out["units"] = units.upper()
    return out


def _csv_columns(path):
    rows = []
    for r in csv.DictReader(open(path)):
        # a run continued with --recover can repeat the last rows before the checkpoint: keep the newest
        try:
            t = float(r.get("time", "nan"))
            while rows and float(rows[-1]["time"]) >= t:
                rows.pop()
        except (TypeError, ValueError):
            pass
        rows.append(r)
    cols = {}
    for k in rows[0].keys():
        try:
            cols[k] = np.array([float(r[k]) if r[k] not in ("", None) else np.nan for r in rows])
        except ValueError:
            pass
    return cols


def read_moose(prefix, field_pressure=True, well=None):
    """A MOOSE run written by em2ex2moose.py + porous_flow-opt, converted to the same quantities as
    read_eclipse(). `prefix` is the -o prefix (it may include a folder)."""
    meta_path = prefix + "_meta.json"
    meta = json.load(open(meta_path)) if os.path.exists(meta_path) else {}
    # surface densities: from the meta file written by em2ex2moose.py (Croft's values if it is missing)
    rho_g = meta.get("rho_gas_sc") or 0.7987446576791414
    rho_w = meta.get("rho_water_sc") or 1020.344079994504
    wells = meta.get("wells") or []
    well = well or (wells[0] if wells else None)

    c = _csv_columns(prefix + "_csv.csv")
    if well is None:   # older runs without a meta file: take the first well found in the CSV
        well = next((k[:-len("_gas_rate_kgps")] for k in c if k.endswith("_gas_rate_kgps")), None)
    t = c["time"]
    ok = c.get("gas_in_place_kg", np.ones_like(t)) != 0      # the first row is written before the postprocessors run
    out = {"t_days": t[ok] / 86400.0, "meta": meta}
    if "gas_in_place_kg" in c:
        out["gip_sm3"] = c["gas_in_place_kg"][ok] / rho_g
    if "water_in_place_kg" in c:
        out["wip_sm3"] = c["water_in_place_kg"][ok] / rho_w
    else:
        # finite-volume inputs report the cumulative aquifer influx instead (the wells produce no water)
        # (minus the water the wells produce, when they produce any)
        cum_in = [k for k in c if k.startswith("aq") and k.endswith("_cum_in_kg")]
        cum_w = [k for k in c if k.endswith("_cum_water_actual_kg")]
        if cum_in:
            out["wip_sm3"] = (sum(c[k] for k in cum_in) + sum(c[k] for k in cum_w))[ok] / rho_w
    if well and "%s_water_rate_actual_kgps" % well in c:   # connection wells producing water by mobility
        out["water_rate_sm3d"] = -c["%s_water_rate_actual_kgps" % well][ok] / rho_w * 86400.0
        out["water_prod_sm3"] = -c["%s_cum_water_actual_kg" % well][ok] / rho_w
    rate_col = "%s_gas_rate_kgps" % well if well else None
    if well and "%s_gas_rate_actual_kgps" % well in c:   # connection wells: the rate actually produced
        rate_col = "%s_gas_rate_actual_kgps" % well
    if rate_col and rate_col in c:
        q = -c[rate_col] / rho_g * 86400.0                       # production positive, sm3/day
        dt = np.diff(t, prepend=t[0])
        out["gas_rate_sm3d"] = q[ok]
        out["gas_cum_sm3"] = np.cumsum(q * dt / 86400.0)[ok]     # rate at the end of each step x step length
    pcol = "%s_pressure" % well if well else None
    gcol = "%s_gas_pressure" % well if well else None
    wcol = "%s_wbp" % well if well else None
    if wcol and wcol in c:          # Eclipse-style WBP (CF-weighted, at the BHP reference depth)
        out["p_well_block_pa"] = np.where(c[wcol][ok] > 0, c[wcol][ok], np.nan)
    elif gcol and gcol in c:        # gas-phase pressure (with capillary pressure), as Eclipse's WBP
        out["p_well_block_pa"] = np.where(c[gcol][ok] > 0, c[gcol][ok], np.nan)   # 0 = not yet computed
    elif pcol and pcol in c:
        out["p_well_block_pa"] = np.where(c[pcol][ok] > 0, c[pcol][ok], np.nan)
    if field_pressure and os.path.exists(prefix + "_exodus.e") and os.path.exists(prefix + ".e"):
        tx, px = hcpv_pressure(prefix + "_exodus.e", prefix + ".e")
        out["t_field_days"], out["p_field_pa"] = tx / 86400.0, px
    return out


def _names(ds, var):
    return ["".join(ch.decode() for ch in row if ch) for row in ds[var][:]]


def hcpv_pressure(exodus_path, mesh_path):
    """Hydrocarbon-pore-volume weighted pressure (as Eclipse FPR) at every MOOSE output time."""
    import netCDF4
    m = netCDF4.Dataset(mesh_path)
    m.set_auto_mask(False)
    nb = m.dimensions["num_el_blk"].size
    conn = np.concatenate([m["connect%d" % b][:].astype(np.int64) - 1 for b in range(1, nb + 1)])
    ev = _names(m, "name_elem_var")
    k = ev.index("pore_volume") + 1
    pv = np.concatenate([m["vals_elem_var%deb%d" % (k, b)][0] for b in range(1, nb + 1)])
    m.close()
    e = netCDF4.Dataset(exodus_path)
    e.set_auto_mask(False)
    if "name_nod_var" not in e.variables or "pwater" not in _names(e, "name_nod_var"):
        # finite-volume output: one value per cell, element variables in the same block order as the mesh
        ev2 = _names(e, "name_elem_var")
        pname = "pgas" if "pgas" in ev2 else "pwater"      # Eclipse's FPR is the gas-phase pressure
        ip, isg = ev2.index(pname) + 1, ev2.index("sgas") + 1
        t = e["time_whole"][:]
        p_out = np.empty(len(t))
        for s in range(len(t)):
            p = np.concatenate([e["vals_elem_var%deb%d" % (ip, b)][s] for b in range(1, nb + 1)])
            sg = np.clip(np.concatenate([e["vals_elem_var%deb%d" % (isg, b)][s] for b in range(1, nb + 1)]), 0, 1)
            w = pv * sg
            p_out[s] = (p * w).sum() / max(w.sum(), 1e-30)
        e.close()
        return t, p_out
    nv = _names(e, "name_nod_var")
    ip, isg = nv.index("pwater") + 1, nv.index("sgas") + 1
    t = e["time_whole"][:]
    p_out = np.empty(len(t))
    for s in range(len(t)):
        p = e["vals_nod_var%d" % ip][s][conn].mean(1)
        sg = np.clip(e["vals_nod_var%d" % isg][s][conn].mean(1), 0, 1)
        w = pv * sg
        p_out[s] = (p * w).sum() / max(w.sum(), 1e-30)
    e.close()
    return t, p_out


# ----------------------------------------------------------------------------------------------
# plotting and tables
# ----------------------------------------------------------------------------------------------
PANELS = [
    # key, title, y label, SI -> plot unit factor
    ("gip_sm3", "Gas in place", "Bscf", 1 / (1e9 * SCF)),
    ("gas_cum_sm3", "Cumulative gas produced", "Bscf", 1 / (1e9 * SCF)),
    ("gas_rate_sm3d", "Gas production rate", "MMscf/day", 1 / (1e6 * SCF)),
    ("p_field_pa", "Field pressure (HCPV-weighted)", "psia", 1 / PSI),
    ("p_well_block_pa", "Well-block pressure", "psia", 1 / PSI),
    ("wip_change", "Water influx into the reservoir", "million stb", 1 / (1e6 * STB)),
    ("water_rate_sm3d", "Well water production rate", "stb/day", 1 / STB),
    ("water_prod_sm3", "Cumulative water produced by the well", "stb", 1 / STB),
]


def _series(d, key):
    if key == "wip_change":
        if "wip_sm3" not in d:
            return None, None
        return d["t_days"], d["wip_sm3"] - d["wip_sm3"][0]
    if key not in d:
        return None, None
    t = d["t_field_days"] if key == "p_field_pa" and "t_field_days" in d else d["t_days"]
    return t, d[key]


def plot_comparison(ecl, runs, panels=None, fig=None, eclipse_style=None, run_styles=None, title=None):
    """Plot Eclipse (black) and each MOOSE run (dashed) on a grid of panels.

    ecl      dict from read_eclipse()
    runs     {label: dict from read_moose()}
    panels   list of keys from PANELS to show (default: all six)
    Returns (fig, axes dict keyed by panel key), so every axis can be edited afterwards."""
    import matplotlib.pyplot as plt
    keys = panels or [p[0] for p in PANELS
                      if _series(ecl, p[0])[0] is not None or any(_series(d, p[0])[0] is not None for d in runs.values())]
    spec = {p[0]: p for p in PANELS}
    ncol = 3 if len(keys) > 4 else 2 if len(keys) > 1 else 1
    nrow = int(np.ceil(len(keys) / ncol))
    if fig is None:
        fig = plt.figure(figsize=(5.4 * ncol, 4.0 * nrow))
    axes = {}
    est = dict(color="k", lw=2.5, label="Eclipse")
    est.update(eclipse_style or {})
    colours = ["tab:red", "tab:blue", "tab:green", "tab:orange", "tab:purple"]
    for n, key in enumerate(keys):
        _, ttl, ylab, fac = spec[key]
        ax = fig.add_subplot(nrow, ncol, n + 1)
        axes[key] = ax
        t, y = _series(ecl, key)
        if t is not None:
            ax.plot(t, y * fac, **est)
        for r, (lab, d) in enumerate(runs.items()):
            t, y = _series(d, key)
            if t is None:
                continue
            st = dict(color=colours[r % len(colours)], lw=2, ls="--", label=lab)
            st.update((run_styles or {}).get(lab, {}))
            if key in ("gas_rate_sm3d", "water_rate_sm3d"):
                st.setdefault("drawstyle", "steps-pre")
            ax.plot(t, y * fac, **st)
        ax.set_title(ttl)
        ax.set_ylabel(ylab)
        ax.set_xlabel("days since start")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig, axes


def _at(d, key, day):
    t, y = _series(d, key)
    if t is None or len(t) == 0:
        return np.nan
    return float(np.interp(day, t, y))


def summary_table(ecl, runs, days=(0, 180, 365, 730, 1000, 1328.5)):
    """Text table of the main quantities at the given days (Eclipse and each run)."""
    cols = [("gip_sm3", "GIP Bscf", 1 / (1e9 * SCF), "%9.3f"), ("gas_cum_sm3", "Gp Bscf", 1 / (1e9 * SCF), "%9.3f"),
            ("p_field_pa", "FPR psia", 1 / PSI, "%9.0f"), ("p_well_block_pa", "WBP psia", 1 / PSI, "%9.0f"),
            ("wip_change", "dWIP MMstb", 1 / (1e6 * STB), "%9.3f"),
            ("water_prod_sm3", "Wp stb", 1 / STB, "%9.1f")]
    width = max([len("Eclipse")] + [len(k) for k in runs]) + 2
    lines = ["%-8s %-*s" % ("day", width, "") + " ".join("%9s" % c[1] for c in cols)]
    for day in days:
        for n, (lab, d) in enumerate([("Eclipse", ecl)] + list(runs.items())):
            vals = " ".join((fmt % (_at(d, k, day) * f)) if not np.isnan(_at(d, k, day)) else "%9s" % "-"
                            for k, _, f, fmt in cols)
            lines.append("%-8s %-*s" % (("%g" % day) if n == 0 else "", width, lab) + vals)
    return "\n".join(lines)


def differences(ecl, run, days=None):
    """Run minus Eclipse for each quantity (plot units) at the given days, as a dict of arrays."""
    days = np.asarray(days if days is not None else ecl["t_days"])
    out = {"day": days}
    for key, _, ylab, fac in PANELS:
        te, ye = _series(ecl, key)
        tm, ym = _series(run, key)
        if te is None or tm is None:
            continue
        out[key + " [" + ylab + "]"] = (np.interp(days, tm, ym) - np.interp(days, te, ye)) * fac
    return out


# ----------------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="Compare MOOSE runs from em2ex2moose.py with the Eclipse results")
    ap.add_argument("--eclipse-h5", required=True, help="Eclipse results (.h5 written by Petrel)")
    ap.add_argument("--run", action="append", required=True, metavar="LABEL=PREFIX",
                    help="MOOSE run: label and the -o prefix given to em2ex2moose.py (repeat for several runs)")
    ap.add_argument("--units", default=None, help="Eclipse unit system (FIELD or METRIC); default: from the run's meta file")
    ap.add_argument("-o", "--output", default="moose_vs_eclipse.png", help="figure file (png, pdf, svg)")
    ap.add_argument("--no-field-pressure", action="store_true", help="skip reading the large MOOSE exodus file")
    ap.add_argument("--days", type=float, nargs="+", default=[0, 180, 365, 730, 1000, 1328.5])
    ap.add_argument("--title", default="MOOSE vs Eclipse")
    a = ap.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    runs = {}
    for r in a.run:
        lab, _, prefix = r.partition("=")
        if not prefix:
            lab, prefix = os.path.basename(r), r
        print("reading %s ..." % prefix)
        runs[lab] = read_moose(prefix, field_pressure=not a.no_field_pressure)
    units = a.units or next((d["meta"].get("units_in") for d in runs.values() if d["meta"].get("units_in")), "FIELD")
    ecl = read_eclipse(a.eclipse_h5, units=units)
    fig, _ = plot_comparison(ecl, runs, title=a.title)
    fig.savefig(a.output, dpi=130)
    print("wrote", a.output)
    print(summary_table(ecl, runs, a.days))


if __name__ == "__main__":
    main()

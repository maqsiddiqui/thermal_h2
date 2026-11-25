"""Split the aquifer influx and pressure of a finite-volume (--fv) MOOSE run between the
separately connected regions of the grid (blocks cut off from each other by faults, where
Eclipse connects them through non-neighbour connections that em2ex does not create).

usage (Miniforge prompt, in the em2ex_to_moose folder):
    python fault_block_split.py runs\\croft_full_fv
    python fault_block_split.py runs\\croft_full_fv --days 0 30 90 180 365 730 1328.5

PREFIX is the -o given to em2ex2moose.py; the script reads PREFIX_exodus.e (results),
PREFIX.i (aquifer coefficient) and PREFIX_csv.csv (MOOSE's own total influx, as a check).
It prints a table and writes PREFIX_fault_blocks.csv and PREFIX_fault_blocks.png.
"""
import argparse
import csv
import re

import numpy as np
import netCDF4
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

PSI = 6894.757
STB_KG = 0.158987 * 1020.34      # kg of water in one stb (only for the printed table)
HEX_FACES = [[0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7], [0, 3, 2, 1], [4, 5, 6, 7]]


def names(d, v):
    return ["".join(c.decode() for c in r if c) for r in d[v][:]]


def regions(C):
    """Label cells by face-connected region (cells sharing a whole face are connected)."""
    faces = np.sort(C[:, HEX_FACES].reshape(-1, 4), axis=1)
    owner = np.repeat(np.arange(len(C)), 6)
    o = np.lexsort(faces.T[::-1])
    fs, ow = faces[o], owner[o]
    same = np.all(fs[1:] == fs[:-1], axis=1)
    a, b = ow[:-1][same], ow[1:][same]
    n = len(C)
    return connected_components(coo_matrix((np.ones(len(a)), (a, b)), shape=(n, n)), directed=False)


def aquifer_coefficient(inp, aq=1):
    """(kg/s/m2/Pa, aquifer time constant in s or None) from the aqN_influx material in the input file."""
    txt = open(inp).read()
    m = re.search(r"\[aq%d_influx\].*?expression\s*=\s*'([^']*)'" % aq, txt, re.S)
    if not m:
        raise SystemExit("could not find the aq%d_influx expression in %s" % (aq, inp))
    expr = m.group(1)
    mult = float(re.match(r"\s*([0-9.eE+-]+)", expr).group(1))
    tau = re.search(r"exp\(-max\(dt, 1\) / ([0-9.eE+-]+)\)", expr)
    return mult, (float(tau.group(1)) if tau else None)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix")
    ap.add_argument("--days", type=float, nargs="*", default=[0, 30, 90, 180, 365, 730, 1000, 1328.5])
    a = ap.parse_args()

    e = netCDF4.Dataset(a.prefix + "_exodus.e")
    e.set_auto_mask(False)
    nb = e.dimensions["num_el_blk"].size
    C = np.vstack([e["connect%d" % (b + 1)][:].astype(np.int64) - 1 for b in range(nb)])
    X = np.c_[e["coordx"][:], e["coordy"][:], e["coordz"][:]]
    t = e["time_whole"][:]
    ev = names(e, "name_elem_var")
    need = ("pwater", "sgas", "p_aq1")
    miss = [v for v in need if v not in ev]
    if miss:
        raise SystemExit("results file lacks %s - is this an --fv run?" % miss)

    def elem(var, step):
        k = ev.index(var) + 1
        return np.concatenate([e["vals_elem_var%deb%d" % (k, b + 1)][step] for b in range(nb)])

    # pore volume: from the mesh written by em2ex2moose.py (same element order)
    m = netCDF4.Dataset(a.prefix + ".e")
    m.set_auto_mask(False)
    mv = names(m, "name_elem_var")
    k = mv.index("pore_volume") + 1
    pv = np.concatenate([m["vals_elem_var%deb%d" % (k, b + 1)][0] for b in range(nb)])
    if len(pv) != len(C):
        raise SystemExit("mesh and results have different numbers of cells")

    nreg, lab = regions(C)
    print("%d separately connected regions" % nreg)

    # aquifer faces (element, area)
    ss = names(e, "ss_names")
    els, area = [], []
    for k, n in enumerate(ss):
        if n.startswith("aquifer_1"):
            el = e["elem_ss%d" % (k + 1)][:] - 1
            sd = e["side_ss%d" % (k + 1)][:] - 1
            for q, s in zip(el, sd):
                P = X[C[q][HEX_FACES[s]]]
                els.append(q)
                area.append(0.5 * np.linalg.norm(np.cross(P[2] - P[0], P[3] - P[1])))
    els, area = np.array(els), np.array(area)
    mult, tau = aquifer_coefficient(a.prefix + ".i")
    reg_aq = lab[els]

    # per region: share of pore volume, initial gas, aquifer face area
    sg0 = np.clip(elem("sgas", 0), 0, 1)
    R = range(nreg)
    big = [r for r in R if pv[lab == r].sum() > 1e-3 * pv.sum()]
    print("\nregion  cells   pore volume  initial gas PV  aquifer faces (area share)")
    for r in big:
        s = lab == r
        print("%6d %7d %12.1f%% %14.1f%% %8d (%.1f%%)" % (
            r, s.sum(), 100 * pv[s].sum() / pv.sum(), 100 * (pv * sg0)[s].sum() / max((pv * sg0).sum(), 1e-30),
            (reg_aq == r).sum(), 100 * area[reg_aq == r].sum() / area.sum()))

    # time series: influx (kg/s, instantaneous, without the step-average factor) and PV-weighted
    # water pressure per region
    # as MOOSE computes it: each step's influx uses that step's pressures and, with a Fetkovich
    # aquifer, the step-average factor tau/dt (1 - exp(-dt/tau)); results must be written every step
    q = np.zeros((len(t), nreg))
    p = np.zeros((len(t), nreg))
    for s_ in range(len(t)):
        pw, pa = elem("pwater", s_), elem("p_aq1", s_)
        dt = max(t[s_] - t[s_ - 1], 1.0) if s_ > 0 else 1.0
        f = tau / dt * (1.0 - np.exp(-dt / tau)) if tau else 1.0
        qf = mult * f * area * (pa[els] - pw[els])
        q[s_] = np.bincount(reg_aq, weights=qf, minlength=nreg)
        p[s_] = np.bincount(lab, weights=pv * pw, minlength=nreg) / np.maximum(np.bincount(lab, weights=pv, minlength=nreg), 1e-30)
    cum = np.vstack([np.zeros(nreg), np.cumsum(q[1:] * np.diff(t)[:, None], axis=0)])

    # check against MOOSE's own total
    try:
        rows = list(csv.reader(open(a.prefix + "_csv.csv")))
        h = rows[0]
        A = np.array(rows[1:], float)
        tc, cm = A[:, h.index("time")], A[:, h.index("aq1_cum_in_kg")]
        moose_total = np.interp(t, tc, cm)
    except Exception:
        moose_total = None

    days = t / 86400
    print("\nday      " + "".join("  region %d: cum influx MMstb  avg p psia" % r for r in big)
          + ("   all regions   MOOSE total" if moose_total is not None else ""))
    for d in a.days:
        if d > days[-1] + 1e-6:
            break
        i = int(np.argmin(abs(days - d)))
        line = "%-8g " % d
        for r in big:
            line += "  %24.3f %11.0f" % (cum[i, r] / STB_KG / 1e6, p[i, r] / PSI)
        if moose_total is not None:
            line += "   %11.3f %11.3f" % (cum[i].sum() / STB_KG / 1e6, moose_total[i] / STB_KG / 1e6)
        print(line)

    with open(a.prefix + "_fault_blocks.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["day"] + ["influx_kgps_r%d" % r for r in big] + ["cum_influx_kg_r%d" % r for r in big]
                   + ["avg_pwater_pa_r%d" % r for r in big])
        for i in range(len(t)):
            w.writerow([days[i]] + list(q[i, big]) + list(cum[i, big]) + list(p[i, big]))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
        for r in big:
            lbl = "region %d (%.0f%% of aquifer faces)" % (r, 100 * area[reg_aq == r].sum() / area.sum())
            ax[0].plot(days, q[:, r], label=lbl)
            ax[1].plot(days, cum[:, r] / STB_KG / 1e6, label=lbl)
            ax[2].plot(days, p[:, r] / PSI, label="region %d" % r)
        ax[1].plot(days, cum.sum(1) / STB_KG / 1e6, "k--", lw=1, label="sum")
        ax[0].set_title("Aquifer influx rate")
        ax[0].set_ylabel("kg/s")
        ax[1].set_title("Cumulative influx")
        ax[1].set_ylabel("million stb")
        ax[2].set_title("Average water pressure (PV-weighted)")
        ax[2].set_ylabel("psia")
        for x in ax:
            x.set_xlabel("days since start")
            x.grid(alpha=0.3)
            x.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(a.prefix + "_fault_blocks.png", dpi=130)
        print("\nwrote %s_fault_blocks.csv and %s_fault_blocks.png" % (a.prefix, a.prefix))
    except Exception as ex:
        print("\nwrote %s_fault_blocks.csv (no plot: %s)" % (a.prefix, ex))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Compare the storage metrics of two (or more) H2 runs, cycle by cycle and in time.

    python compare_h2_runs.py --run "no hysteresis=croft_h2_full" --run "hysteresis=croft_h2_full_hyst" -o compare_hyst

Each --run is "label=prefix": <prefix>_cycles.csv (from postprocess_h2_fv.py) is required,
<prefix>_fv_csv.csv (the MOOSE CSV) is used for the time plots if it exists (or give label=prefix,csvfile).
Writes <o>_cycles.png, <o>_time.png and <o>_table.csv (per-cycle metrics side by side, with differences).
"""
import argparse
import csv
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DAY = 86400.0
PSI = 6894.757293168
COLORS = ["#2a6fdb", "#d9622b", "#2e9e4f", "#8a4fd1"]


def read(fn):
    with open(fn) as f:
        r = list(csv.DictReader(f))
    return {k: np.array([float(x[k]) for x in r]) for k in r[0]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, help="label=prefix[,moose_csv]")
    ap.add_argument("-o", default="compare_h2")
    a = ap.parse_args()
    runs = []
    for s in a.run:
        lab, rest = s.split("=", 1)
        pre, _, ts = rest.partition(",")
        cyc = read(pre + "_cycles.csv")
        ts = ts or (pre + "_fv_csv.csv")
        tsd = read(ts) if os.path.exists(ts) else None
        runs.append((lab, cyc, tsd))

    # ---------------------------------------------------------------- per-cycle panels
    panels = [
        ("h2_recovery_cycle", "H2 recovery per cycle", "fraction of H2 injected", 1),
        ("purity_cycle_mole", "Gas purity (cycle average)", "H2 mole fraction", 1),
        ("h2_mole_frac_end_wd", "Purity at the end of withdrawal", "H2 mole fraction", 1),
        ("working_gas_MMscf", "Working gas volume", "MMscf", 1),
        ("deliverability_start_wd_MMscfd", "Deliverability, start of withdrawal", "MMscf/d", 1),
        ("deliverability_end_wd_MMscfd", "Deliverability, end of withdrawal", "MMscf/d", 1),
        ("bhp_max_inj_psia", "Highest BHP during injection", "psia", 1),
        ("bhp_min_wd_psia", "Lowest BHP during withdrawal", "psia", 1),
        ("h2_in_place_end_t", "H2 left in the reservoir (end of cycle)", "t", 1),
    ]
    fig, axs = plt.subplots(3, 3, figsize=(15, 11))
    w = 0.8 / len(runs)
    for ax, (k, title, unit, _) in zip(axs.flat, panels):
        for i, (lab, cyc, _) in enumerate(runs):
            x = cyc["cycle"] + (i - (len(runs) - 1) / 2) * w
            ax.bar(x, cyc[k], width=w, color=COLORS[i % 4], label=lab)
        lo = min(np.nanmin(c[k]) for _, c, _ in runs)
        hi = max(np.nanmax(c[k]) for _, c, _ in runs)
        if lo > 0 and (hi - lo) < 0.5 * hi:
            ax.set_ylim(lo - 0.15 * (hi - lo) - 1e-9, hi + 0.15 * (hi - lo) + 1e-9)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(unit)
        ax.set_xlabel("cycle")
        ax.set_xticks(runs[0][1]["cycle"])
        ax.grid(alpha=0.3, axis="y")
    axs.flat[0].legend(fontsize=8)
    fig.suptitle("Croft H2 storage: " + " vs ".join(r[0] for r in runs))
    fig.tight_layout()
    fig.savefig(a.o + "_cycles.png", dpi=110)

    # ---------------------------------------------------------------- time plots
    if all(r[2] is not None for r in runs):
        fig, axs = plt.subplots(4, 1, figsize=(11, 13), sharex=True)
        for i, (lab, _, d) in enumerate(runs):
            t = d["time"] / DAY
            c = COLORS[i % 4]
            axs[0].plot(t, d["bhp_Pa"] / PSI, color=c, lw=1, label=lab + ": BHP")
            axs[0].plot(t, d["p_hcpv_Pa"] / PSI, color=c, lw=1, ls="--", label=lab + ": reservoir")
            wd = d["gas_wd_target_kgps"] > 0
            y = np.where(wd, d["h2_molefrac_withdrawn"], np.nan)
            axs[1].plot(t, y, color=c, lw=1.2, label=lab)
            axs[2].plot(t, d["deliverability_sm3pd"] / 28316.85 , color=c, lw=1, label=lab)
            axs[3].plot(t, d["h2_in_place_kg"] / 1e3, color=c, lw=1.2, label=lab)
        axs[0].set_ylabel("psia"); axs[0].set_title("Bottom-hole (solid) and reservoir (dashed) pressure")
        axs[1].set_ylabel("H2 mole fraction"); axs[1].set_title("Purity of the withdrawn gas")
        axs[2].set_ylabel("MMscf/d"); axs[2].set_title("Deliverability (rate possible at the minimum BHP)")
        axs[3].set_ylabel("t"); axs[3].set_title("H2 in place"); axs[3].set_xlabel("days since start of storage")
        for ax in axs:
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(a.o + "_time.png", dpi=110)

    # ---------------------------------------------------------------- table
    keys = [k for k, *_ in panels] + ["h2_withdrawn_t", "ng_withdrawn_t", "aquifer_net_influx_t"]
    base = runs[0][1]
    with open(a.o + "_table.csv", "w", newline="") as f:
        wr = csv.writer(f)
        head = ["cycle"]
        for k in keys:
            head += ["%s [%s]" % (k, lab) for lab, *_ in runs] + ["%s change [%s - %s]" % (k, lab, runs[0][0])
                                                                 for lab, *_ in runs[1:]]
        wr.writerow(head)
        for j, cy in enumerate(base["cycle"]):
            row = [int(cy)]
            for k in keys:
                vals = [r[1][k][j] if j < len(r[1][k]) else np.nan for r in runs]
                row += ["%.6g" % v for v in vals] + ["%.6g" % (v - vals[0]) for v in vals[1:]]
            wr.writerow(row)
    # ---------------------------------------------------------------- console summary
    print("%-34s" % "metric" + "".join("%26s" % ("%s c1 -> c%d" % (lab, len(c["cycle"]))) for lab, c, _ in runs))
    for k, title, unit, _ in panels:
        print("%-34s" % (title[:33]) + "".join("%26s" % ("%.4g -> %.4g" % (c[k][0], c[k][-1])) for _, c, _ in runs))
    print("wrote %s_cycles.png%s, %s_table.csv" % (a.o, ", %s_time.png" % a.o if all(r[2] is not None for r in runs) else "", a.o))


if __name__ == "__main__":
    main()

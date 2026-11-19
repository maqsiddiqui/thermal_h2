#!/usr/bin/env python
"""Per-cycle table, time plots and plume maps of a gen_h2_fv.py run.

    python postprocess_h2_fv.py croft_h2_fv_csv.csv [--exodus croft_h2_fv_exodus.e --mesh croft_h2_fv.e] -o croft_h2
        [--inj-days 182.625 --wd-days 91.3125]

Writes <o>_cycles.csv, <o>_results.png and (with --exodus) <o>_plume.png. Prints the mass-balance check.
"""
import argparse, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PSI = 6894.757293168
DAY = 86400.0
MMSCF = 28316.846592  # sm3
RHO_SC_H2, RHO_SC_NG = 0.08527, 0.7987  # kg/sm3 (H2 at 15 C, 1 atm; native gas from the deck DENSITY)
M_H2, M_NG = 2.016, 18.92


def load(path):
    rows = list(csv.DictReader(open(path)))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def cycles(d, T_INJ, T_WD):
    P = T_INJ + T_WD
    t = d["time"] / DAY
    at = lambda tt, k: np.interp(tt, t, d[k])
    out = []
    for c in range(int(np.floor(t[-1] / P + 1e-6))):
        a, b, e = c * P, c * P + T_INJ, (c + 1) * P
        inj = (t > a) & (t <= b + 1e-6)
        wd = (t > b) & (t <= e + 1e-6)
        h2i = at(b, "cum_h2_injected_kg") - at(a, "cum_h2_injected_kg")
        h2w = at(e, "cum_h2_withdrawn_kg") - at(b, "cum_h2_withdrawn_kg")
        ngw = at(e, "cum_ng_withdrawn_kg") - at(b, "cum_ng_withdrawn_kg")
        aq = at(e, "aq1_cum_in_kg") - at(a, "aq1_cum_in_kg")
        y = d["h2_molefrac_withdrawn"][wd]
        wdr = d["withdrawal_sm3pd"][wd] if "withdrawal_sm3pd" in d else np.array([])
        dpot = d["deliverability_sm3pd"][wd] if "deliverability_sm3pd" in d else np.array([])
        wgv_h2, wgv_ng = h2w / RHO_SC_H2, ngw / RHO_SC_NG
        out.append(dict(
            cycle=c + 1,
            p_well_start_psia=at(a, "p_well_Pa") / PSI,
            p_well_end_inj_psia=at(b, "p_well_Pa") / PSI,
            p_well_max_psia=d["p_well_Pa"][inj | wd].max() / PSI,
            p_well_end_wd_psia=at(e, "p_well_Pa") / PSI,
            p_res_start_psia=at(a, "p_hcpv_Pa") / PSI,
            p_res_end_inj_psia=at(b, "p_hcpv_Pa") / PSI,
            p_res_end_wd_psia=at(e, "p_hcpv_Pa") / PSI,
            h2_injected_t=h2i / 1e3, h2_withdrawn_t=h2w / 1e3, ng_withdrawn_t=ngw / 1e3,
            h2_recovery_cycle=h2w / h2i if h2i > 0 else np.nan,
            h2_recovery_cumulative=at(e, "cum_h2_withdrawn_kg") / max(at(e, "cum_h2_injected_kg"), 1e-30),
            h2_mass_frac_withdrawn=h2w / max(h2w + ngw, 1e-30),
            h2_mole_frac_start_wd=y[0] if len(y) else np.nan,
            h2_mole_frac_end_wd=y[-1] if len(y) else np.nan,
            h2_in_place_end_t=at(e, "h2_in_place_kg") / 1e3,
            ng_in_place_end_t=at(e, "ng_in_place_kg") / 1e3,
            aquifer_net_influx_t=aq / 1e3,
            working_gas_MMscf=(wgv_h2 + wgv_ng) / MMSCF,
            working_gas_h2_MMscf=wgv_h2 / MMSCF,
            working_gas_ng_MMscf=wgv_ng / MMSCF,
            h2_injected_MMscf=h2i / RHO_SC_H2 / MMSCF,
            purity_cycle_mole=(h2w / M_H2) / max(h2w / M_H2 + ngw / M_NG, 1e-30),
            deliverability_start_wd_MMscfd=(dpot[0] / MMSCF) if len(dpot) else np.nan,
            deliverability_end_wd_MMscfd=(dpot[-1] / MMSCF) if len(dpot) else np.nan,
            withdrawal_rate_end_MMscfd=(wdr[-1] / MMSCF) if len(wdr) else np.nan,
            injectivity_end_inj_MMscfd=(np.interp(b, t, d["injectivity_sm3pd"]) / MMSCF) if "injectivity_sm3pd" in d else np.nan,
            bhp_max_inj_psia=(d["bhp_Pa"][inj].max() / PSI) if "bhp_Pa" in d and inj.any() else np.nan,
            bhp_min_wd_psia=(d["bhp_Pa"][wd].min() / PSI) if "bhp_Pa" in d and wd.any() else np.nan,
        ))
    return out


def metrics(d, cyc, a, t, P):
    """the three evaluation metrics: deliverability, working gas volume, purity"""
    fig, ax = plt.subplots(3, 1, figsize=(9, 11))
    wd = d["gas_wd_target_kgps"] > 0 if "gas_wd_target_kgps" in d else d["gas_wd_rate_kgps"] > 0
    inj = d["h2_inj_target_kgps"] > 0 if "h2_inj_target_kgps" in d else d["h2_inj_rate_kgps"] > 0
    for a_ in ax[[0, 2]]:
        for c in range(int(np.ceil(t[-1] / P))):
            a_.axvspan(c * P, c * P + a.inj_days, color="tab:blue", alpha=0.06, lw=0)
            a_.axvspan(c * P + a.inj_days, (c + 1) * P, color="tab:orange", alpha=0.08, lw=0)
        a_.grid(alpha=0.3)
    if "deliverability_sm3pd" in d:
        ax[0].plot(t, d["deliverability_sm3pd"] / MMSCF, color="tab:green", lw=1.4,
                   label="deliverability: rate possible at the minimum BHP")
        if "injectivity_sm3pd" in d:
            ax[0].plot(t, d["injectivity_sm3pd"] / MMSCF, color="tab:blue", ls="--", lw=1.2,
                       label="injectivity: H2 rate possible at the maximum BHP (initial pressure)")
        ax[0].plot(t, np.where(wd, d["withdrawal_sm3pd"] / MMSCF, np.nan), color="k", lw=1.6, label="actual withdrawal rate")
        ax[0].plot(t, np.where(inj, d["h2_inj_rate_kgps"] * DAY / RHO_SC_H2 / MMSCF, np.nan), color="tab:blue", lw=1.2,
                   label="actual H2 injection rate")
        ax[0].set_ylabel("MMscf/d"); ax[0].set_yscale("log"); ax[0].legend(fontsize=8)
    ax[0].set_title("Deliverability (blue band = injection, orange = withdrawal)")
    if cyc:
        n = np.array([c["cycle"] for c in cyc])
        h = np.array([c["working_gas_h2_MMscf"] for c in cyc]); g = np.array([c["working_gas_ng_MMscf"] for c in cyc])
        hi = np.array([c["h2_injected_MMscf"] for c in cyc])
        ax[1].bar(n, h, 0.6, color="tab:green", label="H2 withdrawn")
        ax[1].bar(n, g, 0.6, bottom=h, color="tab:gray", label="native gas withdrawn")
        ax[1].plot(n, hi, "o", color="tab:blue", label="H2 injected")
        for k in range(len(n)):
            ax[1].text(n[k], h[k] + g[k], "%.3g" % (h[k] + g[k]), ha="center", va="bottom", fontsize=8)
        ax[1].set_xticks(n); ax[1].set_xlabel("cycle"); ax[1].set_ylabel("MMscf"); ax[1].legend(fontsize=8)
        ax[1].grid(alpha=0.3, axis="y")
    ax[1].set_title("Working gas volume per cycle (gas withdrawn, standard conditions)")
    ax[2].plot(t, np.where(wd, d["h2_molefrac_withdrawn"], np.nan), color="tab:green", lw=1.4, label="withdrawn gas, instantaneous")
    if cyc:
        tm = np.array([(c["cycle"] - 1) * P + a.inj_days + a.wd_days / 2 for c in cyc])
        ax[2].plot(tm, [c["purity_cycle_mole"] for c in cyc], "s", color="k", label="cycle average")
    ax[2].set_ylim(0, 1.02); ax[2].set_ylabel("H2 mole fraction"); ax[2].set_xlabel("days since start of storage")
    ax[2].legend(fontsize=8); ax[2].set_title("Gas purity")
    fig.tight_layout(); fig.savefig(a.o + "_metrics.png", dpi=110)
    print("wrote %s_metrics.png" % a.o)


def plume(exo, mesh, out, title):
    import netCDF4
    from scipy.spatial import cKDTree
    nm = lambda ds, v: ["".join(c.decode() for c in r if c) for r in ds[v][:]]

    def geo(ds):
        X = np.c_[ds["coordx"][:], ds["coordy"][:], ds["coordz"][:]]
        nb = ds.dimensions["num_el_blk"].size
        C = np.vstack([ds["connect%d" % (b + 1)][:] - 1 for b in range(nb)])
        return X, C, nb

    def ev(ds, name, step, nb):
        k = nm(ds, "name_elem_var").index(name)
        return np.concatenate([ds["vals_elem_var%deb%d" % (k + 1, b + 1)][step] for b in range(nb)])
    m = netCDF4.Dataset(mesh); e = netCDF4.Dataset(exo)
    Xm, Cm, nbm = geo(m); Xe, Ce, nbe = geo(e)
    idx = cKDTree(Xe[Ce].mean(1)).query(Xm[Cm].mean(1))[1]   # mesh element -> exodus element
    I = ev(m, "cell_i", 0, nbm).astype(int); J = ev(m, "cell_j", 0, nbm).astype(int)
    col = I * 10000 + J
    ucol, inv = np.unique(col, return_inverse=True)
    first = np.array([np.nonzero(inv == c)[0][0] for c in range(len(ucol))])
    polys = [Xm[Cm[i][:4], :2] for i in first]
    t = e["time_whole"][:] / DAY
    picks = sorted(set([1, 2, len(t) - 2, len(t) - 1])) if len(t) > 2 else list(range(len(t)))
    fig, ax = plt.subplots(2, len(picks), figsize=(4.6 * len(picks), 8), squeeze=False)
    sg_any = np.zeros(len(ucol))
    for s in picks:
        np.maximum.at(sg_any, inv, ev(e, "sgas", s, nbe)[idx])
    cx = np.array([p.mean(0) for p in polys])[sg_any > 1e-3]
    lim = (cx[:, 0].min() - 250, cx[:, 0].max() + 250, cx[:, 1].min() - 250, cx[:, 1].max() + 250) if len(cx) else None
    from matplotlib.collections import PolyCollection
    for c, s in enumerate(picks):
        sg = ev(e, "sgas", s, nbe)[idx]; xh = ev(e, "xh2", s, nbe)[idx]
        v1 = np.full(len(ucol), np.nan); v2 = np.zeros(len(ucol))
        np.fmax.at(v1, inv, np.where(sg > 1e-3, xh, np.nan))
        np.maximum.at(v2, inv, sg)
        for r, (v, lab, cmap, cl) in enumerate(((v1, "max H2 mass fraction of gas in column", "viridis", (0, 1)),
                                                (v2, "max gas saturation in column", "magma_r", (0, 0.9)))):
            pc = PolyCollection(polys, array=np.ma.masked_invalid(v), cmap=cmap, edgecolors="none")
            pc.set_clim(*cl)
            a = ax[r, c]; a.add_collection(pc); a.autoscale(); a.set_aspect("equal")
            if lim:
                a.set_xlim(lim[0], lim[1]); a.set_ylim(lim[2], lim[3])
            a.set_title("%s\nday %.0f" % (lab, t[s]), fontsize=8)
            plt.colorbar(pc, ax=a, shrink=0.6)
    fig.suptitle(title); fig.tight_layout(); fig.savefig(out, dpi=100)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv"); ap.add_argument("-o", default="croft_h2")
    ap.add_argument("--exodus"); ap.add_argument("--mesh")
    ap.add_argument("--inj-days", type=float, default=182.625); ap.add_argument("--wd-days", type=float, default=91.3125)
    ap.add_argument("--title", default="Croft H2 storage (FV)")
    a = ap.parse_args()
    d = load(a.csv); t = d["time"] / DAY; P = a.inj_days + a.wd_days
    cyc = cycles(d, a.inj_days, a.wd_days)
    with open(a.o + "_cycles.csv", "w") as f:
        if cyc:
            f.write(",".join(cyc[0]) + "\n")
            for c in cyc:
                f.write(",".join("%.6g" % c[k] for k in c) + "\n")
    fig, ax = plt.subplots(4, 1, figsize=(9, 13), sharex=True)
    for a_ in ax:
        for c in range(int(np.ceil(t[-1] / P))):
            a_.axvspan(c * P, c * P + a.inj_days, color="tab:blue", alpha=0.06, lw=0)
            a_.axvspan(c * P + a.inj_days, (c + 1) * P, color="tab:orange", alpha=0.08, lw=0)
        a_.grid(alpha=0.3)
    if "bhp_Pa" in d:
        b = np.where((d["h2_inj_target_kgps"] > 0) | (d["gas_wd_target_kgps"] > 0), d["bhp_Pa"], np.nan)
        ax[0].plot(t, b / PSI, "-", color="tab:purple", lw=1.4, label="bottom-hole pressure (reference depth)")
        ax[0].axhline(d["bhp_max_limit_Pa"][-1] / PSI, color="tab:purple", ls=":", lw=1, label="injection limit = initial reservoir pressure")
        ax[0].axhline(d["bhp_min_limit_Pa"][-1] / PSI, color="tab:purple", ls="-.", lw=1, label="withdrawal limit (minimum BHP)")
    ax[0].plot(t, d["p_well_Pa"] / PSI, "k-", label="gas pressure, top completion cell")
    ax[0].plot(t, d["p_hcpv_Pa"] / PSI, "r-", label="reservoir (gas-PV weighted)")
    ax[0].plot(t, d["aq1_pressure_at_datum"] / PSI, "b--", label="aquifer at datum")
    ax[0].set_ylabel("psia"); ax[0].legend(fontsize=8); ax[0].set_title("Pressure (blue = injection, orange = withdrawal)")
    ax[1].plot(t, d["h2_in_place_kg"] / 1e3, "g-", label="H2 in place")
    ax[1].plot(t, d["cum_h2_injected_kg"] / 1e3, "g:", label="cumulative H2 injected")
    ax[1].plot(t, d["cum_h2_withdrawn_kg"] / 1e3, "g--", label="cumulative H2 withdrawn")
    ax[1].plot(t, d["cum_ng_withdrawn_kg"] / 1e3, "m--", label="cumulative native gas withdrawn")
    ax[1].set_ylabel("tonnes"); ax[1].legend(fontsize=8); ax[1].set_title("Inventories")
    wd = d["gas_wd_rate_kgps"] > 0
    ax[2].plot(t, np.where(wd, d["h2_molefrac_withdrawn"], np.nan), "g-", lw=1.6)
    ax[2].set_ylim(0, 1.02); ax[2].set_ylabel("H2 mole fraction"); ax[2].set_title("Purity of the withdrawn gas")
    ax[3].plot(t, d["aq1_cum_in_kg"] / 1e3, "b-")
    ax[3].set_ylabel("tonnes"); ax[3].set_xlabel("days since start of storage")
    ax[3].set_title("Cumulative aquifer influx (+ in, - out)")
    fig.tight_layout(); fig.savefig(a.o + "_results.png", dpi=110)
    metrics(d, cyc, a, t, P)
    print("wrote %s_cycles.csv, %s_results.png" % (a.o, a.o))
    for c in cyc:
        print("cycle %2d: well %.0f -> %.0f -> %.0f psia, reservoir %.0f -> %.0f -> %.0f psia; H2 in %.1f t, out %.1f t "
              "(recovery %.2f), native gas out %.1f t; withdrawn H2 mole fraction %.3f -> %.3f; working gas %.1f MMscf, "
              "cycle purity %.3f, deliverability %.2f -> %.2f MMscf/d"
              % (c["cycle"], c["p_well_start_psia"], c["p_well_end_inj_psia"], c["p_well_end_wd_psia"], c["p_res_start_psia"],
                 c["p_res_end_inj_psia"], c["p_res_end_wd_psia"], c["h2_injected_t"], c["h2_withdrawn_t"], c["h2_recovery_cycle"],
                 c["ng_withdrawn_t"], c["h2_mole_frac_start_wd"], c["h2_mole_frac_end_wd"], c["working_gas_MMscf"],
                 c["purity_cycle_mole"], c["deliverability_start_wd_MMscfd"], c["deliverability_end_wd_MMscfd"]))
    i = -1
    print("mass balance at day %.2f: H2 %.3g kg (of %.4g injected), native gas %.3g kg (of %.4g), water %.3g kg (of %.4g)"
          % (t[i], d["h2_balance_kg"][i], d["cum_h2_injected_kg"][i], d["ng_balance_kg"][i], d["ng_in_place_kg"][0],
             d["water_balance_kg"][i], d["water_in_place_kg"][0]))
    print("xng range over the run: %.4f - %.4f" % (d["min_xng"].min(), d["max_xng"].max()))
    if a.exodus and a.mesh:
        plume(a.exodus, a.mesh, a.o + "_plume.png", a.title)


if __name__ == "__main__":
    main()

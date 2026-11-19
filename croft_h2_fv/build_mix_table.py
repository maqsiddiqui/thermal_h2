"""Gas-mixture property table for an H2 storage model: rho(p, w_ng), mu(p, w_ng) at reservoir T.
w_ng = mass fraction of native (cushion) gas in the gas phase, stored in the 'temperature' column
because MOOSE's TabulatedBicubicFluidProperties interpolates in (p, T).
Native gas: the deck's PVDG table (above its top pressure: DAK/Sutton, scaled to PVDG).
H2: CoolProp. Density: Amagat additive volumes, times a non-ideal-mixing correction
rho_GERG(CH4-H2)/rho_Amagat(CH4-H2) at the same H2 mole fraction.
Viscosity: Herning-Zipperer mixing of native-gas and H2 viscosities.
    python build_mix_table.py DECK.DATA out.csv [T_kelvin]     (T fitted to PVDG if omitted)"""
import sys, numpy as np, CoolProp.CoolProp as CP
from pvt_tools import load_pvdg, gas_mw, fit_temperature, native_gas_functions
DECK = sys.argv[1]; OUT = sys.argv[2] if len(sys.argv) > 2 else "h2_mix_pvt.csv"
p, bg, mu, rho_sc = load_pvdg(DECK)
if len(sys.argv) > 3:
    T = float(sys.argv[3]); print("T = %.2f K (given)" % T)
else:
    T, e = fit_temperature(p, bg, rho_sc); print("T = %.2f K (%.1f C) fitted to PVDG, rms %.2f %%" % (T, T - 273.15, 100 * e))
M_NG = gas_mw(rho_sc); M_H2 = 2.016e-3; M_C1 = 16.043e-3
print("native gas: rho_sc %.4f kg/m3, M %.2f g/mol, PVDG %.2f-%.2f MPa" % (rho_sc, M_NG * 1e3, p.min() / 1e6, p.max() / 1e6))
rho_ng, mu_ng = native_gas_functions(p, bg, mu, rho_sc, T)
P = np.unique(np.r_[np.linspace(1.5e6, 20e6, 38), np.linspace(21e6, 35e6, 15)])
W = np.r_[np.linspace(0, 0.7, 15), np.linspace(0.72, 0.9, 10), np.linspace(0.91, 1, 10)]
rows = []
AS = CP.AbstractState("HEOS", "Methane&Hydrogen"); AS.specify_phase(CP.iphase_gas)
for Pi in P:
    rh = CP.PropsSI("D", "P", Pi, "T", T, "Hydrogen"); mh = CP.PropsSI("V", "P", Pi, "T", T, "Hydrogen")
    rc = CP.PropsSI("D", "P", Pi, "T", T, "Methane")
    rn, mn = rho_ng(Pi)[0], mu_ng(Pi)[0]
    for w in W:
        yng = (w / M_NG) / (w / M_NG + (1 - w) / M_H2)            # native-gas mole fraction
        rho_a = 1 / (w / rn + (1 - w) / rh)
        corr = 1.0
        if 0 < w < 1:  # same H2 mole fraction in the CH4-H2 analogue
            wc = yng * M_C1 / (yng * M_C1 + (1 - yng) * M_H2)
            AS.set_mole_fractions([yng, 1 - yng]); AS.update(CP.PT_INPUTS, Pi, T); rg = AS.rhomass()
            corr = rg * (wc / rc + (1 - wc) / rh)
        sq_ng, sq_h = np.sqrt(M_NG), np.sqrt(M_H2)
        mu_m = (yng * mn * sq_ng + (1 - yng) * mh * sq_h) / (yng * sq_ng + (1 - yng) * sq_h)
        rows.append((Pi, w, rho_a * corr, mu_m))
# extend linearly to w = -0.05 and 1.05 so Newton iterates slightly outside [0, 1] see smooth properties
ext = []
for Pi in P:
    blk = [r for r in rows if r[0] == Pi]
    (w0, r0, m0), (w1, r1, m1) = blk[0][1:], blk[1][1:]
    (wa, ra, ma), (wb, rb, mb) = blk[-2][1:], blk[-1][1:]
    ext += [(Pi, -0.05, r0 - 0.05 * (r1 - r0) / (w1 - w0), m0 - 0.05 * (m1 - m0) / (w1 - w0))] + blk + \
           [(Pi, 1.05, rb + 0.05 * (rb - ra) / (wb - wa), mb + 0.05 * (mb - ma) / (wb - wa))]
with open(OUT, "w") as f:
    f.write("pressure,temperature,density,viscosity\n")
    for r in ext: f.write("%.6e,%.6f,%.8e,%.8e\n" % r)
print("wrote", OUT, len(P), "pressures x", len(W) + 2, "fractions")

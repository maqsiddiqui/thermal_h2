"""PVT helpers for the H2 storage models (native gas from the Eclipse PVDG table)."""
import numpy as np
import sys
import os
# eclio.py comes with em2ex_to_moose: set EM2EX_TO_MOOSE to that folder, or keep it next to this folder
_here = os.path.dirname(os.path.abspath(__file__))
for _p in (os.environ.get("EM2EX_TO_MOOSE", ""), os.path.join(_here, "..", "em2ex_to_moose"), _here):
    if _p and os.path.exists(os.path.join(_p, "eclio.py")):
        sys.path.insert(0, _p)
        break
import eclio

PSI = 6894.757293168
R = 8.314462618
T_SC = 288.706       # 60 F
P_SC = 101352.9      # 14.7 psia
M_AIR = 28.9647e-3


def z_dak(ppr, tpr):
    """Dranchuk & Abou-Kassem (1975) Z-factor, Newton on reduced density."""
    A = [0.3265, -1.0700, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210]
    z = 1.0
    for _ in range(200):
        rr = 0.27 * ppr / (z * tpr)
        f = (1 + (A[0] + A[1] / tpr + A[2] / tpr**3 + A[3] / tpr**4 + A[4] / tpr**5) * rr
             + (A[5] + A[6] / tpr + A[7] / tpr**2) * rr**2
             - A[8] * (A[6] / tpr + A[7] / tpr**2) * rr**5
             + A[9] * (1 + A[10] * rr**2) * (rr**2 / tpr**3) * np.exp(-A[10] * rr**2))
        zn = 0.5 * z + 0.5 * f
        if abs(zn - z) < 1e-10:
            break
        z = zn
    return z


def sutton_pc(sg):
    """Sutton (1985) pseudo-criticals for natural gas: returns (Ppc [Pa], Tpc [K])."""
    ppc = (756.8 - 131.07 * sg - 3.6 * sg**2) * PSI
    tpc = (169.2 + 349.5 * sg - 74.0 * sg**2) / 1.8
    return ppc, tpc


def load_pvdg(deck_path):
    d = eclio.Deck(deck_path)
    p, bg, mu = d.pvdg()
    rho_sc = d.densities()[2]
    return p, bg, mu, rho_sc


def gas_mw(rho_sc):
    """Molar mass [kg/mol] of the native gas from its standard-condition density."""
    rho_air = P_SC * M_AIR / (R * T_SC)
    return rho_sc / rho_air * M_AIR


def rho_dak(p, T, rho_sc):
    """Native-gas density [kg/m3] from DAK + Sutton for gas gravity rho_sc/rho_air."""
    M = gas_mw(rho_sc); sg = M / M_AIR
    ppc, tpc = sutton_pc(sg)
    return np.array([pp * M / (z_dak(pp / ppc, T / tpc) * R * T) for pp in np.atleast_1d(p)])


def fit_temperature(p, bg, rho_sc, pmin=2e6):
    """Reservoir temperature [K] that makes DAK/Sutton reproduce the PVDG Bg curve (PVDG has no T)."""
    sel = p >= pmin
    target = rho_sc / bg[sel]
    Ts = np.arange(280.0, 420.0, 0.25)
    err = [np.sqrt(np.mean((rho_dak(p[sel], T, rho_sc) / target - 1) ** 2)) for T in Ts]
    i = int(np.argmin(err))
    return Ts[i], err[i]


def native_gas_functions(p, bg, mu, rho_sc, T):
    """rho(P), mu(P) of the native gas: PVDG inside its range; beyond the top PVDG pressure, DAK density
    scaled to match PVDG at the last point and viscosity extrapolated linearly."""
    lp = np.log(p); lr = np.log(rho_sc / bg)
    pmax = p.max(); scale = (rho_sc / bg[-1]) / rho_dak(pmax, T, rho_sc)[0]
    slope_mu = (mu[-1] - mu[-2]) / (p[-1] - p[-2])
    def rho(P):
        P = np.atleast_1d(P).astype(float)
        out = np.exp(np.interp(np.log(P), lp, lr))
        hi = P > pmax
        if hi.any():
            out[hi] = scale * rho_dak(P[hi], T, rho_sc)
        return out
    def visc(P):
        P = np.atleast_1d(P).astype(float)
        out = np.interp(P, p, mu)
        hi = P > pmax
        out[hi] = mu[-1] + slope_mu * (P[hi] - pmax)
        return out
    return rho, visc

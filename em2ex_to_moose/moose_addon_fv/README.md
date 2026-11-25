# Add-on for `--fv` inputs: two-point flux, Eclipse rock compressibility, tabulated rel-perm with hysteresis, fault connections

Five objects, ten files, all needed by inputs written with `em2ex2moose.py --fv` (the third one unless
`--relperm corey` is given, the fourth and fifth unless `--nnc off` is given).

## 1. FVPorousFlowAdvectiveFluxTPFA – two-point flux for cell-centred PorousFlow

Needed by inputs written with `em2ex2moose.py --fv`. It replaces MOOSE's `FVPorousFlowAdvectiveFlux`,
which on Croft's dipping corner-point grid is not in hydrostatic balance (a reservoir at rest starts to
flow) and upwinds permeability; with it the 19-layer Croft model diverges in the first step.

On each internal face:

    flux out of cell E = sum over phases  (X kr rho / mu)_upwind * T * Phi
    T   = 1 / (1/tE + 1/tN),   t = k_dir (n . D) / (D . D)     Eclipse's corner-point transmissibility
    Phi = (pE - pN) + rho_face g . (xN - xE)   gravity along the same line as the pressure difference

D runs from the cell centre to the face centre; k_dir is the permeability along the grid direction of the
connection (vertical permeability for mostly-horizontal faces, otherwise the horizontal permeability along
D). On the full Croft grid this reproduces Eclipse's TRANX, TRANY and TRANZ (ratio 1.000 over thousands of
sampled connections; taking n.K.n instead overestimates the vertical ones 2.2x on the dipping layers). A
hydrostatic state gives exactly zero flux on any grid.

## 2. ADPorousFlowPorosityRock – Eclipse's ROCK compressibility

    porosity = phi_ref (1 + X + X^2/2),   X = c_r (P - p_ref)

with P the effective fluid pressure, phi_ref the deck porosity at the ROCK reference pressure p_ref and
c_r the ROCK compressibility (Eclipse's own expression). It also works on MOOSE versions older than
the AD version of PorousFlowPorosity (e.g. April 2026).

## 3. ADPorousFlowRelativePermeabilityHysteresisTable – Eclipse rel-perm tables and hysteresis

Relative permeability of one phase from the deck's SWFN/SGFN table (linear interpolation, constant beyond
the ends, as Eclipse), with Eclipse's hysteresis options (SATOPTS HYSTER, EHYSTR item 2):

    hysteresis = none        kr = m kr_drainage(S)
    hysteresis = imbibition  kr = m kr_imbibition(S)            water with EHYSTR 1 or 3: in every cell
    hysteresis = killough    gas with EHYSTR 2 or 3 (Killough 1976, as ECLIPSE and OPM Flow):
        S >= S_hy  kr = kr_d(S)                                  on the drainage curve
        S <  S_hy  S_t = S_dc + (S_hy - S_dc) / (1 + C (S_hy - S_dc))      trapped gas (Land)
                   C   = 1/(S_ic - S_dc) - 1/(S_max - S_dc)
                   S_n = S_ic + (S - S_t)(S_max - S_ic)/(S_hy - S_t)
                   kr  = kr_i(S_n) kr_d(S_hy) / kr_i(S_max)       0 below S_t

S_dc, S_ic are the critical gas saturations of the drainage (SATNUM) and imbibition (IMBNUM) tables, S_max
the end of the drainage table, S_hy the largest gas saturation the cell has had: an aux variable `sg_max`
updated with `max(sg_max, sgas)` at the end of every time step (lagged, like MOOSE's own hysteresis
materials, so the curve choice never switches inside a Newton iteration). m (`kr_multiplier`) carries
Eclipse's end-point scaling of the maximum water rel-perm (KRW, or IKRW for the imbibition curve).
Checked against a Python implementation of OPM Flow's formula along a drainage–imbibition–redrainage path
(difference 5e-15).

## 4. FVPorousFlowNNCFlux + NNCRelationshipManager – Eclipse's connections across faults

Where a fault has throw, the cells on its two sides do not share a mesh face, so MOOSE has no face to put a
flux on and the fault is sealed. Eclipse connects them (TRANX/TRANY over the overlapping part, plus its NNCs).
On Croft a fault right next to CROFT_1 (between I = 36 and 37) was sealed this way, which made the well block
pressure 36-160 psi too low in the hysteresis run.

The converter writes `<prefix>_nnc.csv` with every Eclipse connection (TRANX/Y/Z from the INIT, NNC1/NNC2
from the EGRID with TRANNNC) whose two MOOSE cells share no mesh face, Eclipse's transmissibility in SI
(m3; merged layers summed). `FVPorousFlowNNCFlux` (one per component equation) adds to cell a

    1/V_a  sum_b  T_ab  sum_p  (X kr rho / mu)_p,upwind * [ (p + o_p)_a - (p + o_p)_b + rho_p,avg g.(x_b - x_a) ]

with p the primary pressure (implicit at both cells) and o_p a per-phase offset (0 for water, Pc for gas);
kr, rho, mu come from auxiliary variables of the start of the step (as the well connections), the mass
fractions may be implicit functors. Centroids are the FV ones (vertex average), so a hydrostatic state gives
no flow. `NNCRelationshipManager` (added automatically) ghosts the partner cells and adds the Jacobian
entries, so the connections work in parallel (checked: identical residuals and Newton history on 1, 2 and 4
processes). The input sets `allow_renumbering = false` in `[Mesh]`: the file refers to element ids in the
Exodus order.

## Add them to your app (once, on every machine that runs these inputs, including the HPC)

    cd ~/projects/thermal_h2
    mkdir -p include/fvkernels src/fvkernels
    mkdir -p include/materials src/materials
    cp FVPorousFlowAdvectiveFluxTPFA.h include/fvkernels/
    cp FVPorousFlowAdvectiveFluxTPFA.C src/fvkernels/
    cp ADPorousFlowPorosityRock.h include/materials/
    cp ADPorousFlowPorosityRock.C src/materials/
    cp ADPorousFlowRelativePermeabilityHysteresisTable.h include/materials/
    cp ADPorousFlowRelativePermeabilityHysteresisTable.C src/materials/
    mkdir -p include/relationshipmanagers src/relationshipmanagers
    cp FVPorousFlowNNCFlux.h include/fvkernels/
    cp FVPorousFlowNNCFlux.C src/fvkernels/
    cp NNCRelationshipManager.h include/relationshipmanagers/
    cp NNCRelationshipManager.C src/relationshipmanagers/
    make -j 4
    ./thermal_h2-opt --registry | grep -E "FVPorousFlowAdvectiveFluxTPFA|ADPorousFlowPorosityRock|HysteresisTable|NNC"   # five lines

(Your app must include the PorousFlow module, as it does for the hydrogen add-on.) Then run the FV inputs
with your app instead of `porous_flow-opt`:

    mpiexec -n 12 ~/projects/thermal_h2/thermal_h2-opt -i croft_full_fv.i

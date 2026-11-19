# C++ add-ons for the finite-volume Croft models

| files | object | used by | what it does |
|---|---|---|---|
| `FVPorousFlowAdvectiveFluxTPFA.h/.C` | `FVPorousFlowAdvectiveFluxTPFA` | history and H2 FV inputs | two-point flux with Eclipse's corner-point transmissibility, hydrostatically consistent (unchanged from `em2ex_to_moose/moose_addon_fv`) |
| `ADPorousFlowPorosityRock.h/.C` | `ADPorousFlowPorosityRock` | history and H2 FV inputs | Eclipse ROCK compressibility (unchanged) |
| `PorousFlowMassFractionLinear.h/.C` | `PorousFlowMassFractionLinear`, `ADPorousFlowMassFractionLinear` | H2 inputs | mass fraction = multiplier × variable + offset (the 1e-5 dissolved-gas fraction in water) |
| `FVMaxIncrementDamper.h/.C` | `FVMaxIncrementDamper` | H2 FV input | Newton damping for finite-volume variables: largest change per iteration, optional bounds, optional skipping of cells where another variable (gas saturation) is small; `verbose = true` prints the limiting cell (run MOOSE with `--keep-cout` to see all ranks) |

Install into the PorousFlow module (as on the HPC):

    cd moose/modules/porous_flow
    mkdir -p include/dampers src/dampers
    cp FVPorousFlowAdvectiveFluxTPFA.h include/fvkernels/;  cp FVPorousFlowAdvectiveFluxTPFA.C src/fvkernels/
    cp ADPorousFlowPorosityRock.h include/materials/;       cp ADPorousFlowPorosityRock.C src/materials/
    cp PorousFlowMassFractionLinear.h include/materials/;   cp PorousFlowMassFractionLinear.C src/materials/
    cp FVMaxIncrementDamper.h include/dampers/;             cp FVMaxIncrementDamper.C src/dampers/
    make -j 8

or into your own app (same sub-folders of the app's `include/` and `src/`). All objects register under
`PorousFlowApp`, so they are available in any app that includes the PorousFlow module.

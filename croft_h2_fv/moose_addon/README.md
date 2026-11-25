# C++ add-ons for the finite-volume Croft models

| files | object | used by | what it does |
|---|---|---|---|
| `FVPorousFlowAdvectiveFluxTPFA.h/.C` | `FVPorousFlowAdvectiveFluxTPFA` | history and H2 FV inputs | two-point flux with Eclipse's corner-point transmissibility, hydrostatically consistent (unchanged from `em2ex_to_moose/moose_addon_fv`) |
| `ADPorousFlowPorosityRock.h/.C` | `ADPorousFlowPorosityRock` | history and H2 FV inputs | Eclipse ROCK compressibility (unchanged) |
| `PorousFlowMassFractionLinear.h/.C` | `PorousFlowMassFractionLinear`, `ADPorousFlowMassFractionLinear` | H2 inputs | mass fraction = multiplier × variable + offset (the 1e-5 dissolved-gas fraction in water) |
| `ADPorousFlowRelativePermeabilityHysteresisTable.h/.C` | `ADPorousFlowRelativePermeabilityHysteresisTable` | history and H2 FV inputs made with `--relperm deck` (default) | tabulated rel-perm (Eclipse SWFN/SGFN, linear) with Eclipse's hysteresis: Killough scanning curves for gas between the SATNUM and IMBNUM tables, imbibition table for water (EHYSTR option 3); per-cell multiplier for IKRW; historical maximum Sg from a lagged aux variable `sg_max` |
| `FVMaxIncrementDamper.h/.C` | `FVMaxIncrementDamper` | H2 FV input | Newton damping for finite-volume variables: largest change per iteration, optional bounds, optional skipping of cells where another variable (gas saturation) is small; `verbose = true` prints the limiting cell (run MOOSE with `--keep-cout` to see all ranks) |
| `FVPorousFlowNNCFlux.h/.C`, `NNCRelationshipManager.h/.C` | `FVPorousFlowNNCFlux`, `NNCRelationshipManager` | history and H2 FV inputs whose history was made with fault connections (`--nnc on`, default) | flow between cells across faults with throw, which share no mesh face, using Eclipse's transmissibilities from `<prefix>_nnc.csv`; the relationship manager ghosts the partner cells so it works in parallel |

Install into the PorousFlow module (as on the HPC):

    cd moose/modules/porous_flow
    mkdir -p include/dampers src/dampers
    cp FVPorousFlowAdvectiveFluxTPFA.h include/fvkernels/;  cp FVPorousFlowAdvectiveFluxTPFA.C src/fvkernels/
    cp ADPorousFlowPorosityRock.h include/materials/;       cp ADPorousFlowPorosityRock.C src/materials/
    cp PorousFlowMassFractionLinear.h include/materials/;   cp PorousFlowMassFractionLinear.C src/materials/
    cp FVMaxIncrementDamper.h include/dampers/;             cp FVMaxIncrementDamper.C src/dampers/
    cp ADPorousFlowRelativePermeabilityHysteresisTable.h include/materials/
    cp ADPorousFlowRelativePermeabilityHysteresisTable.C src/materials/
    mkdir -p include/relationshipmanagers src/relationshipmanagers
    cp FVPorousFlowNNCFlux.h include/fvkernels/;              cp FVPorousFlowNNCFlux.C src/fvkernels/
    cp NNCRelationshipManager.h include/relationshipmanagers/; cp NNCRelationshipManager.C src/relationshipmanagers/
    make -j 8

or into your own app (same sub-folders of the app's `include/` and `src/`). All objects register under
`PorousFlowApp`, so they are available in any app that includes the PorousFlow module.

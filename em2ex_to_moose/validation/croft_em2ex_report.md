# em2ex2moose conversion report

Source: `/tmp/claude-0/-home-claude/f8745664-a354-5c1f-b624-a5c17e7f48a5/scratchpad/croft_in/CROFT_files/HM_MODEL.DATA`

- Layers merged exactly as in the em2ex mesh (--coarsen-z 5: 92 -> 19 layers). Restart saturations, SWL and region arrays not in the mesh are pore-volume weighted over the same groups; wells and aquifer connections remapped.
- Porosity passed to MOOSE is em2ex's PORO x NTG (Eclipse pore-volume convention).
- Horizontal permeability multiplied by NTG (Eclipse applies NTG to X/Y transmissibility).
- Cells with effective porosity <= 0.0001 are excluded: 0 active cells.
- Mesh geometry from em2ex: croft_19.e (same nodes and elements; z converted to elevation, x/y shifted by the origin below).
- Porosity, NTG, PERMX/Y/Z and SATNUM from the em2ex mesh's element variables.
- Mesh: 48805 HEX8 elements, 56026 nodes, 0 elements with mirrored node order (left-handed grid).
- 2090 of 141534 logical internal faces are non-conforming (fault throw); MOOSE treats them as no-flow.
- MAPAXES present (2142150.2 18811116.0 2142150.2 18814398.0 2145431.0 18814398.0); coordinates are used as stored in COORD (no transform applied).
- Coordinates shifted by origin (0.0, 0.0) m; z = -depth.
- Pore volume in MOOSE mesh: 4.73454e+07 m3 = 100.0001% of Eclipse PORV (4.73454e+07 m3).
- Water PVT from PVTW: Bw=1.0320 at 1.972e+07 Pa, cw=4.39e-10 1/Pa, mu=0.000283 Pa s -> SimpleFluidProperties (reservoir density 988.7 kg/m3).
- Gas PVT from PVDG (14 points, 2.76e+06-2.07e+07 Pa) tabulated to croft_em2ex_gas_pvt.csv for TabulatedBicubicFluidProperties (isothermal: same values at every T).
- 6 SATNUM regions in use -> one MOOSE block per SATNUM (satnum_1, satnum_2, satnum_3, satnum_4, satnum_5, satnum_6), each with its own Corey fit.
- Rel-perm satnum_1 <- Corey fit to SWFN/SGFN table 1: Swr=0.720 nw=3.93 krw_max=0.903, Sgr=0.000 ng=3.00 krg_max=0.331.
- Rel-perm satnum_2 <- Corey fit to SWFN/SGFN table 2: Swr=0.726 nw=4.00 krw_max=1.000, Sgr=0.000 ng=3.00 krg_max=0.325.
- Rel-perm satnum_3 <- Corey fit to SWFN/SGFN table 3: Swr=0.471 nw=4.00 krw_max=1.000, Sgr=0.000 ng=3.00 krg_max=0.574.
- Rel-perm satnum_4 <- Corey fit to SWFN/SGFN table 4: Swr=0.266 nw=1.49 krw_max=0.992, Sgr=0.000 ng=1.00 krg_max=0.762.
- Rel-perm satnum_5 <- Corey fit to SWFN/SGFN table 5: Swr=0.150 nw=1.50 krw_max=1.000, Sgr=0.000 ng=1.00 krg_max=0.865.
- Rel-perm satnum_6 <- Corey fit to SWFN/SGFN table 6: Swr=0.100 nw=1.50 krw_max=1.000, Sgr=0.000 ng=1.00 krg_max=0.909.
- NOTE: the deck uses end-point scaling (SWL varies 0.10-0.73 cell by cell); MOOSE Corey curves use one set of end points per block. Initial saturations still honour the cell values (from the restart).
- Initial saturations from restart HM_MODEL.X0000; initial pressure hydrostatic from EQUIL (MOOSE has Pc = 0, so the restart's capillary-laden phase pressure would not be in equilibrium). Difference to the restart pressure: median -0.21 bar, 99th pct |.| 0.80 bar. Use --restart-pressure to take it as is.
- Initial gas in place in the MOOSE model: 2.12544e+08 sm3 = 7.50593e+06 Mscf (surface conditions, from PVDG).
- Pc = 0 surrogate: 1664 gas-zone cells hold water above the table residual (Eclipse keeps it there with capillary pressure). They are put in 109 sub-blocks (initial Sw binned by 0.025) whose water curve starts at that Sw, so this water is immobile until it is exceeded; gas curves are unchanged. Use --mobile-connate-water to switch this off.
- Aquifer 1: initial pressure defaulted -> equilibrium with connected cells 1.982e+07 Pa.
- Aquifer 1 (FETKOVICH): boundaries aquifer_1_Ip aquifer_1_Jm; flux = 6.09e-10 kg/s/m2/Pa x (p - p_aq(z,t)); Fetkovich PI 1000 (deck units) spread over 4.47e+05 m2; finite: p_aq drops by cumulative influx / (rho V ct), V ct = 0.2306 m3/Pa.
- Wells: CROFT_1. Rates from WCONHIST/WCONPROD/WCONINJE/WCONINJH converted to mass rates with surface densities gas 0.799, water 1020.3 kg/m3; split over connections by Kh/CF.
- Well rates averaged over 30-day periods (cumulative volumes preserved at period ends).
  NOTE: wells are rate-controlled point sources for the controlled phase (gas for GRAT); BHP limits are NOT enforced and water is only produced where a water rate is given.

# Croft hydrogen storage: finite-volume MOOSE model (full grid, HPC)

Cyclic hydrogen storage in Croft, started from the end of the finite-volume history match. It is the
finite-volume (cell-centred, two-point flux) version of the earlier nodal 19-layer storage model, with the
same schedule and physics:

* 10 cycles of: inject **1000 t of pure H2 over 6 months** (182.625 d, 0.06338 kg/s), then withdraw
  **1000 t of gas over 3 months** (91.3125 d, 0.12675 kg/s), with whatever H2/native-gas mixture is in
  the completion cells;
* well CROFT_1 at the crest, cell (45, 43), completions as in the history (4 cells on the full grid,
  rate split by KH);
* start state = pressure and gas saturation of every cell at the end of the history run (1 Jan 2006);
* grid, rock, rel-perm, capillary pressure, rock compressibility, water PVT and the Fetkovich aquifer are
  copied from the history input, so the storage model is exactly the history model;
* the aquifer is re-based to its depleted state at the end of history (Eclipse-style, cell pressures and
  step averaging, as in the history model).

## Files

| file | what it is |
|---|---|
| `gen_h2_fv.py` | writes the storage input from a finite-volume history input and its Exodus output |
| `h2_mix_pvt.csv` | density and viscosity of the native-gas + H2 mixture, 1.5–35 MPa, at 76.4 °C |
| `build_mix_table.py`, `pvt_tools.py` | rebuild that table from the deck (needs CoolProp; reproduces the file exactly) |
| `postprocess_h2_fv.py` | per-cycle table, time plots and plume maps; prints the mass-balance check |
| `run_h2_hpc.sh` | Slurm script |
| `moose_addon/` | the four C++ add-ons the input needs (two new, two from the history model) |
| `checks/` | the 19-layer check runs described below (CSV, figures, tables) |
| `README.md` | this file |

## Physics (what the input does)

* **3 components, 2 phases.** Component 0 = water (equation on `pwater`), 1 = native gas (on `xng`),
  2 = hydrogen (on `sgas`). Phase 0 = water, phase 1 = gas (native gas + H2 mixture). Primary variables,
  one value per cell: `pwater`, `sgas`, `xng` (native-gas mass fraction of the gas).
* **Mixture properties.** ρ(p, xng) and μ(p, xng) from `h2_mix_pvt.csv`: native gas from the deck's PVDG
  (molar mass 18.92 g/mol, temperature 76.4 °C fitted to PVDG), H2 from CoolProp, density by additive
  volumes with a GERG non-ideal-mixing correction, viscosity by Herning–Zipperer. At xng = 1 the table
  equals the history model's gas table (within 0.01 % at 12–20 MPa). MOOSE's tabulated fluid interpolates
  in (p, T), so xng is passed in the temperature slot; the model is isothermal.
* **Dissolved gas (numerical).** The water phase carries 1e-5 of dissolved gas with the local composition.
  This gives xng a conservation equation where there is no free gas, so Newton stays well posed at the
  gas front. It is about 0.6 % of the native gas in place and is counted in all inventories.
* **Mixing.** Advection of the mixture (composition-dependent density and viscosity) plus the numerical
  dispersion of the upwind two-point scheme (about half a cell: ~20 m laterally, ~0.5 m vertically on the
  full grid). No explicit diffusion/dispersion kernels: molecular diffusion is ~3 m in 6 months.
* **Well (rate targets with bottom-hole-pressure limits).** Injection: pure H2 at the target rate, unless
  the bottom-hole pressure (BHP) would exceed **the initial reservoir pressure** (default 2865.9 psia at the
  BHP reference depth 6509.19 ft from WELSPECS: EQUIL 2875 psia at 6653.5 ft moved up the initial gas
  gradient; Eclipse's WBP at day 0 is 2865 psia). Withdrawal: gas at the target rate, unless the BHP would
  fall below a minimum (default: the lowest well-block pressure of the history run, 2147 psia on 19 layers).
  When a limit is reached the rate drops so that the BHP stays exactly at the limit.
  Each completion cell c is a well connection with Eclipse's connection factor CF_c (COMPDAT; the four
  CF are proportional to KH, sum 3.4739 cP·rb/d/psi = 9.27e-13 m³):
  q_c = min(share_c × target, CF_c ρ λ |p_limit(z_c) − p_gas,c|), with λ = krg/μg + krw/μw (injection,
  ρ = H2 density at the cell pressure) or krg/μg (withdrawal, ρ = cell gas density); p_limit(z) is the
  limit carried from the reference depth along the initial gas gradient. Mobilities are taken from the last
  converged step, pressures and composition are implicit. Each connection is limited on its own (a
  connection below its share does not pass the shortfall to the others). Withdrawn gas is split between
  native gas and H2 by each cell's composition, so purity and recovery are results. No wellbore flow
  model and no water production.
* **Aquifer.** Fetkovich, as in the history (Eclipse's step-averaged influx on cell pressures), re-based:
  p0 = p0_history − (cumulative influx at the end of history)/capacity. Water leaving the reservoir carries
  the dissolved gas of the cell it leaves; water entering carries native gas.
* **Same as the history model:** capillary pressure (Brooks–Corey per sub-block, scaled from PCW/SWLPC),
  rel-perm (Corey fits, per block), rock compressibility, permeability, TPFA flux kernel. Gas–water
  capillary pressure and rel-perm are assumed the same for the H2 mixture as for the native gas.

## Putting it in your thermal_h2 repository

On your PC, unzip `croft_h2_fv.zip` inside your local clone of thermal_h2 so that the folder is
`thermal_h2/croft_h2_fv/`, then:

    cd thermal_h2
    git checkout <your branch>            # e.g. the branch that holds the TPFA kernel
    git add croft_h2_fv
    git commit -m "Croft H2 storage: finite-volume generator, add-ons, checks"
    git push

On the HPC: `cd <thermal_h2 clone> && git pull`. The C++ files are then in `croft_h2_fv/moose_addon/`
(copy them into the PorousFlow module, step 1 below) and the scripts in `croft_h2_fv/`.
The check outputs (`checks/`, ~1 MB) can be left out of the repository if you prefer.

## Step-by-step: generating and running the full-grid model on the HPC

You need the full-grid history run you already made on the HPC (`croft_full_fv.i`, its mesh
`croft_full_fv.e` and its output `croft_full_fv_exodus.e` with the last time step at 114,782,400 s).

1. **Add the C++ add-ons to MOOSE (once).** The input needs `FVPorousFlowAdvectiveFluxTPFA` and
   `ADPorousFlowPorosityRock` (you already have both, unchanged) and two new objects:
   `ADPorousFlowMassFractionLinear` (mass fractions with the small dissolved-gas fraction) and
   `FVMaxIncrementDamper` (Newton damping for finite-volume variables):

        source ~/set_env.sh
        cd $SCRATCH3DIR/moose/modules/porous_flow
        mkdir -p include/dampers src/dampers
        cp <thermal_h2>/croft_h2_fv/moose_addon/PorousFlowMassFractionLinear.h include/materials/
        cp <thermal_h2>/croft_h2_fv/moose_addon/PorousFlowMassFractionLinear.C src/materials/
        cp <thermal_h2>/croft_h2_fv/moose_addon/FVMaxIncrementDamper.h include/dampers/
        cp <thermal_h2>/croft_h2_fv/moose_addon/FVMaxIncrementDamper.C src/dampers/
        make -j 8
        ./porous_flow-opt --registry | grep -E "TPFA|PorosityRock|MassFractionLinear|FVMaxIncrementDamper"   # 5 lines

   (or copy them into your own app's `include/` and `src/` folders and `make` it.)

2. **Generate the storage input** next to the history run (Python with numpy, scipy, netCDF4):

        cd <folder with croft_full_fv.i, croft_full_fv.e, croft_full_fv_exodus.e>
        python <thermal_h2>/croft_h2_fv/gen_h2_fv.py --hist-input croft_full_fv.i \
               --hist-exodus croft_full_fv_exodus.e --deck CROFT_files/HM_MODEL.DATA -o croft_h2_fv

   `--deck` gives the connection factors, BHP reference depth and initial pressure (eclio.py from
   em2ex_to_moose is found through `EM2EX_TO_MOOSE=<folder>` or a sibling folder `em2ex_to_moose`;
   without `--deck` the Croft values are used). `--p-max-psia` / `--p-min-psia` override the limits.

   It writes `croft_h2_fv.i`, `croft_h2_fv.e`, copies `h2_mix_pvt.csv` and prints a report
   (`croft_h2_fv_report.md`). Check in the report: the last time step is t = 114,782,400 s; all 236,316
   cells matched; the aquifer re-base (~314 psi lower); 139 capillary-pressure objects copied.

3. **Optional quick check (recommended, ~1 h):** 2 short cycles before the 10-cycle run:

        python .../gen_h2_fv.py --hist-input croft_full_fv.i --hist-exodus croft_full_fv_exodus.e --deck CROFT_files/HM_MODEL.DATA \
               -o croft_h2_test --mass-t 20 --inj-days 2 --wd-days 1 --cycles 2 --dtmax-days 0.25

   Run it like step 4 and check with `postprocess_h2_fv.py croft_h2_test_csv.csv -o test --inj-days 2 --wd-days 1`:
   no failed steps, 20 t withdrawn per cycle, H2 balance within a few kg.

4. **Run** with `sbatch run_h2_hpc.sh` (edit the paths). Expect roughly 3–5 h per cycle on 16 cores,
   i.e. 30–50 h for 10 cycles (scaled from the 19-layer check; the first cycle is the slowest while the
   plume forms). Memory: the direct solver (MUMPS) needs
   about 8 GB for its factors plus ~6 GB for the model, spread over the ranks; ask for the whole node
   (`--mem=100GB`, 12–24 ranks). If MUMPS stops with error −9, regenerate with `--mumps-relax 200`;
   if memory runs out, use `--mumps-ooc` (out-of-core, slower). Checkpoints are written every step;
   continue an interrupted run by adding `--recover` to the command.

5. **Post-process** (on the HPC or after copying the CSV, Exodus and mesh):

        python postprocess_h2_fv.py croft_h2_fv_csv.csv --exodus croft_h2_fv_exodus.e --mesh croft_h2_fv.e -o croft_h2

   It writes `croft_h2_cycles.csv` (per cycle: pressures and BHP, H2 in/out, recovery, the three
   metrics below, inventories, aquifer), `croft_h2_results.png` (pressures with the BHP and its limits,
   inventories, purity, aquifer), `croft_h2_metrics.png` and `croft_h2_plume.png`, and prints the
   mass-balance check.

## Evaluation metrics (`<o>_metrics.png`, `<o>_cycles.csv`)

* **Deliverability**: the gas rate the well can deliver with its BHP at the minimum (standard conditions,
  MMscf/d), at every time step, with the actual withdrawal and injection rates and the injectivity (the
  H2 rate the well could take with its BHP at the initial reservoir pressure). Per cycle: at the start and
  at the end of withdrawal (`deliverability_start_wd_MMscfd`, `deliverability_end_wd_MMscfd`).
* **Working gas volume**: gas withdrawn in each cycle at standard conditions (MMscf), split into H2 and
  native gas, against the H2 injected (`working_gas_MMscf`, `working_gas_h2_MMscf`, `working_gas_ng_MMscf`,
  `h2_injected_MMscf`). Standard densities: H2 0.08527 kg/sm³ (15 °C, 1 atm), native gas 0.7987 kg/sm³
  (deck DENSITY).
* **Gas purity**: H2 mole fraction of the withdrawn gas, instantaneous and averaged over each withdrawal
  (`purity_cycle_mole`).

With the default 1000 t per cycle the BHP stays below the injection limit in the first cycle (it peaks
around 2630 psia against 2866 psia). To find the well's storage capacity, raise `--mass-t` (or shorten
`--inj-days`): injection then runs at the limit and the injected mass becomes a result.

## Outputs (CSV, every time step)

| column | meaning |
|---|---|
| `h2_inj_target_kgps`, `gas_wd_target_kgps` | target rates |
| `h2_inj_rate_kgps`, `gas_wd_rate_kgps` | actual rates (lower than the target when a BHP limit is reached) |
| `bhp_Pa`, `bhp_max_limit_Pa`, `bhp_min_limit_Pa` | bottom-hole pressure at the reference depth and its two limits |
| `deliverability_sm3pd` | gas rate the well could deliver now with the BHP at its minimum (standard m³/d) |
| `injectivity_sm3pd` | H2 rate the well could take now with the BHP at the injection limit (standard m³/d) |
| `withdrawal_sm3pd` | actual withdrawal rate (standard m³/d) |
| `h2_wd_rate_kgps`, `ng_wd_rate_kgps` | H2 and native gas actually withdrawn (sum = gas rate) |
| `cum_h2_injected_kg`, `cum_h2_withdrawn_kg`, `cum_ng_withdrawn_kg` | cumulatives |
| `h2_recovery_cum` | cumulative H2 withdrawn / injected |
| `h2_molefrac_withdrawn` | H2 mole fraction of the withdrawn gas |
| `xng_well`, `p_well_Pa`, `p_well_max_Pa` | composition and gas pressure in the top completion cell |
| `p_hcpv_Pa`, `hcpv_m3` | gas-pore-volume-weighted gas pressure (as Eclipse FPR), gas pore volume |
| `aq1_influx_kgps`, `aq1_cum_in_kg`, `aq1_pressure_at_datum` | aquifer (+ = into the reservoir) |
| `h2_in_place_kg`, `ng_in_place_kg`, `water_in_place_kg` | inventories integrated over the 3D fields |
| `h2_balance_kg`, `ng_balance_kg`, `water_balance_kg` | in place − (initial + in − out): should stay ~0 |
| `min_xng`, `max_xng`, `nl_its` | composition bounds, Newton iterations |

Exodus output (3D) at the start and at the end of every injection and withdrawal period: `pwater`, `pgas`,
`sgas`, `xng`, `xh2` (H2 mass fraction of the gas), `yh2` (H2 mole fraction), `gas_density`, `phi`.

## Checks done here

The full grid needs more memory than the 7.6 GB machine used here, so the physics was checked on the
19-layer grid (same generator, same objects), and the full-grid input was checked for set-up only.

**History base.** The 19-layer FV history (`em2ex2moose.py --fv --capillary deck`, current version with
`ADPorousFlowPorosityRock`) was rerun to its end (1,328.5 d, 90 steps, no failed steps) and used as start.

**1. Switching test** (`checks/croft_fv19_switch_test_csv.csv`, before the pressure limits were added):
2 cycles of 20 t H2 in over 2 days and 20 t gas out over 1 day. Both switches converge, no failed steps,
exactly 20.00 t withdrawn per cycle, H2 balance −0.06 kg of 40 t.

**2. Pressure limits** (2-day/1-day cycles with 20 t targets and limits set so that they bind):
* injection limit lowered to 2500 psia (`checks/croft_fv19_injection_limit_test_*`): from the first step
  the BHP sits at 2500.0 psia and the H2 rate drops from 0.021 to 0.008 kg/s (target 0.116 kg/s);
* minimum BHP raised to 2500 psia (`checks/croft_fv19_withdrawal_limit_test_*`): after half a day of
  withdrawal the BHP reaches 2500.0 psia and stays there while the rate falls from 0.232 to 0.100 kg/s;
  the actual rate then equals the deliverability exactly, as it should;
* mass balances stay within 14 kg (H2) and 114 kg (native gas, of 76,790 t).

**3. First full cycle with the default limits** (`checks/croft_fv19_cycle1_*`; 1000 t in over 182.6 d,
1000 t out over 91.3 d; injection limit 2865.9 psia, minimum BHP 2146.8 psia):

| | FV 19 layers (this model) | nodal 19 layers (earlier model) |
|---|---|---|
| BHP peak during injection (limit 2865.9 psia) | 2,632 psia: limit not reached | (no well model) |
| well-block pressure: start → end of injection → end of withdrawal | 2,486 → 2,617 → 2,452 psia | 2,480 → 2,661 → 2,439 psia (water pressure) |
| reservoir pressure (gas-PV weighted) | 2,552 → 2,623 → 2,543 psia | 2,535 → 2,626 → 2,538 psia |
| H2 withdrawn / native gas withdrawn | 566.6 t / 433.5 t (H2 recovery 57 %) | 399.6 t / 600.4 t (40 %) |
| **working gas volume** (gas withdrawn) | 253.8 MMscf (234.7 H2 + 19.2 native); H2 injected 414.2 MMscf | |
| **purity** (H2 mole fraction of withdrawn gas) | 1.00 → 0.81, cycle average 0.925 | 0.998 → 0.73 |
| **deliverability** at the minimum BHP | 97.5 MMscf/d at the start of withdrawal, 46.8 at the end (actual rate 4.6 → 1.8 MMscf/d) | |
| water pushed into the aquifer over the cycle | 37,100 t | 49,000 t |
| H2 / native-gas mass balance at the end | −0.13 t of 1000 t / +1.7 t of 76,790 t | 0.3 t |

The well is far from both limits with 1000 t per cycle: the target rates control the whole cycle, and the
deliverability (what the well could produce down to 2147 psia) is 20–27 times the 3-month withdrawal
rate. The FV model recovers more H2 at higher purity than the nodal model: its two-point scheme keeps the
H2 plume compact at the crest (one composition per cell, connections through faces only), while the nodal
scheme smears the plume over the nodes it shares with surrounding cells. The composition variable leaves
[0, 1] only in cells with almost no gas (gas saturation below 1e-3, where it is not damped and has no
effect on the gas); it returns within a few steps. Run time: 85 time steps (12 failed and repeated with
smaller steps), 1.7 h on 2 cores.

**4. Full grid (set-up only).** The generator was run on the 92-layer history input with a surrogate end
state (the 19-layer end state mapped onto the fine cells). The input passes `--check-input` and the first
residual is evaluated on all 236,316 cells; the first Jacobian needs more than 7.6 GB. MUMPS estimates
its factors at 7.7 GB (in-core) for the 709,000 unknowns, so on the HPC ask for the whole node.

**Development notes (why the numerical settings are what they are).**
* Without Newton damping, cells at the edge of the H2 plume with almost no gas received Newton updates
  of the composition of several units, the composition went negative (−4), the mixture table clamped,
  and Newton stalled. MOOSE's `MaxIncrement` damper only works for finite-element variables, hence the
  add-on `FVMaxIncrementDamper`: pressure ≤ 10 bar, saturation ≤ 0.2 and composition ≤ 0.2 per
  iteration, composition kept in [0, 1]. Composition is not damped in cells with gas saturation below
  1e-3, where it is poorly determined; damping it there throttled the whole update to a few per cent and
  stalled the run at day 52.
* `nl_rel_tol = 1e-3`, `nl_abs_tol = 1e-4` kg/s: tighter relative tolerances made Newton stall at the
  plume front (linear convergence from upwind switching); the mass balances above show the chosen values
  keep the error below 0.1 %. (The earlier nodal model used 2e-3.)

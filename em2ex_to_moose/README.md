# em2ex → MOOSE: Croft primary depletion and history match

Three steps, all from the Miniforge prompt:

1. **em2ex** makes the mesh from the Eclipse grid.
2. **`em2ex2moose.py`** writes the MOOSE input file (and everything it reads) for that mesh.
3. **`compare_with_eclipse.py`** compares the MOOSE run with Eclipse, from the command line or in the
   notebook `compare_with_eclipse.ipynb`.

| file | what it is |
|---|---|
| `em2ex2moose.py` | the input-file generator (run this) |
| `petrel2moose.py`, `eclio.py`, `cornerpoint.py`, `gsg.py` | the engine it uses to read the Eclipse deck (PVT, rel-perm, aquifer, wells, restart) |
| `compare_with_eclipse.py` | comparison with Eclipse: command line, or functions for the notebook |
| `compare_with_eclipse.ipynb` | the same comparison in Jupyter, with examples of changing the plots |

## One-time setup

```
conda install -c conda-forge numpy scipy netCDF4 h5py matplotlib jupyter
```

## 1. Make the mesh with em2ex

In your em2ex folder (master, which has `--coarsen-z`):

```
python em2ex.py -f --convert-to-m --coarsen-z 5 --coarsen-z-pinch -o croft_19.e croft.grdecl
```

* `--convert-to-m` is required: Croft is in feet, and the generator needs metres.
* `--coarsen-z 5` merges 92 layers into 19, like the validated MOOSE model (48,805 cells). Leave it out for the
  full 92-layer grid (236,316 cells; much slower to run).
* `--coarsen-z-pinch` because the Croft deck uses `PINCH` (for 5-layer groups it changes nothing).
* Don't use `--translate`, `--mapaxes`, `--refine-xy` or `--extract-*`: they move or cut the grid, so it no
  longer lines up with the Eclipse case.

## 2. Write the MOOSE input

Put the em2ex mesh and this folder's scripts wherever you like; give the paths:

```
python em2ex2moose.py --mesh croft_19.e --case CROFT_files\HM_MODEL.DATA --coarsen-z 5 -o croft_em2ex
```

* `--coarsen-z` must be the same value you gave em2ex (omit it if you didn't coarsen).
* `--case` is the Eclipse deck; `HM_MODEL.EGRID`, `.INIT` and `.X0000` must be next to it.
* `-o croft_em2ex` is the prefix of everything written (it can include a folder, e.g. `-o runs\croft_em2ex`).

Takes about 20 s and writes:

| file | contents |
|---|---|
| `croft_em2ex.i` | MOOSE input: two-phase water/gas, PVDG gas, 6 SATNUM Corey sets, Fetkovich aquifer, CROFT_1 rate history (30-day averages) |
| `croft_em2ex.e` | the em2ex mesh – same nodes and elements – plus the initial state, block names and aquifer sidesets the input needs |
| `croft_em2ex_gas_pvt.csv` | gas PVT table read by the input |
| `croft_em2ex_wells.csv` | well connection(s) |
| `croft_em2ex_report.md` | every step and approximation, with the checks (pore volume, gas in place, aquifer, rel-perm fits) |
| `croft_em2ex_meta.json` | unit system, surface densities, wells – read by the comparison script |

What comes from where:

* **em2ex mesh:** geometry, porosity, NTG, PERMX/Y/Z, SATNUM (which also sets the blocks).
* **Eclipse case:** gas/water PVT, rel-perm tables, initial saturations (restart `.X0000`, merged over the
  same layer groups), hydrostatic initial pressure from EQUIL, aquifer (AQUFETP/AQUANCON), well
  connections and rate history.

Useful options: `--older-moose` (if your MOOSE build stops with *unused parameter … time_integration_scheme*), `--average-rates 0` (exact rate schedule instead of 30-day averages), `--restart-step -1`
(initial state from EQUIL only), `--mobile-connate-water`, `--no-checkpoint`. `python em2ex2moose.py -h` lists them.

## 3. Run MOOSE

From the folder that holds the files above (the input refers to them by name):

```
mpiexec -n 2 porous_flow-opt -i croft_em2ex.i
```

Files made on Windows can be copied to Linux or an HPC as they are: the generator writes Linux (LF) line
endings, which MOOSE requires. (Files from an older copy of the generator may have Windows line endings;
fix them with `sed -i 's/\r$//' croft_em2ex.i croft_em2ex_gas_pvt.csv croft_em2ex_wells.csv`.)

About 70 min on 2 cores for the 1,328.5-day history. It writes `croft_em2ex_csv.csv` (time series) and
`croft_em2ex_exodus.e` (3D results, ~470 MB). Checkpoints are on, so a stopped run continues with
`porous_flow-opt -i croft_em2ex.i --recover`.

## 4. Compare with Eclipse

**Command line:**

```
python compare_with_eclipse.py --eclipse-h5 CROFT_files\HM_MODEL.h5 --run "MOOSE (em2ex mesh)=croft_em2ex" -o croft_vs_eclipse.png
```

* `--run LABEL=PREFIX` – PREFIX is the `-o` you gave `em2ex2moose.py`. Repeat `--run` to overlay several runs.
* `--no-field-pressure` skips reading the big `_exodus.e` file (saves a minute or two).
* `--days 0 365 730` sets the days in the printed table.

It saves a 6-panel figure (gas in place, cumulative gas, gas rate, field pressure, well-block pressure,
water influx) and prints a table at chosen days.

**Notebook:** start Jupyter in this folder (`jupyter notebook` or `jupyter lab`), open
`compare_with_eclipse.ipynb`, set the two paths in the first code cell and run all cells. It uses the same
functions, then shows how to pick panels, change styles and limits, plot your own figure from the data,
and print differences from Eclipse.

## Checks (Croft, `--coarsen-z 5`)

* The em2ex elements line up with the Eclipse cells to within 3.4 mm.
* Pore volume 4.73454e7 m³ = Eclipse PORV; initial gas in place 7.506 Bscf (Eclipse FGIP 7.505).
* The generated input file is the same as the validated Croft model's (only the aquifer coefficient differs
  in the 9th digit), and the mesh properties agree with it within 0.15 %.
* Results against Eclipse (full 1,328.5-day history, 4,422 s on 2 cores) are in `validation/`:

| day | | GIP Bscf | Gp Bscf | field p psia | well-block p psia | water influx MMstb |
|---|---|---|---|---|---|---|
| 365 | Eclipse | 3.760 | 3.745 | 2335 | 2244 | 2.385 |
| | MOOSE (em2ex mesh) | 3.813 | 3.747 | 2295 | 2216 | 2.429 |
| 1328.5 | Eclipse | 3.382 | 4.123 | 2522 | 2510 | 3.091 |
| | MOOSE (em2ex mesh) | 3.437 | 4.123 | 2514 | 2502 | 3.203 |

  MOOSE's gas-in-place output reads 0.7 % above Eclipse FGIP from day 0 (as in the validated model), and
  stays about 0.055 Bscf above it throughout. `croft_em2ex_and_validated_vs_eclipse.png` overlays the earlier
  validated MOOSE run: the input files are the same and only the mesh differs. The em2ex-mesh run follows
  Eclipse's well-block pressure more closely (about 30 psi low at the end of production vs about 190 psi low).

## Memory: how fine a grid fits

The input uses MUMPS, a direct linear solver: robust, but its memory grows quickly with the number of
cells. Measured on Croft (2 processes):

| em2ex `--coarsen-z` | layers | cells | MUMPS memory | runs in 7 GB? |
|---|---|---|---|---|
| 5 | 19 | 48,805 | ~1.5 GB | yes (about 75 min) |
| 3 | 31 | 79,629 | ~2.9 GB | yes (a few hours) |
| 2 | 46 | 118,158 | ~5.2 GB | no (killed) |
| none | 92 | 236,316 | ~13 GB | no |

If MUMPS doesn't have enough memory, the log shows `FACTOR_OUTMEMORY` / `DIVERGED_PC_FAILED` and the
residual never drops, or the process is killed. The iterative alternative in the input (`active = asm`)
uses little memory but does not converge on Croft, so for finer grids use a machine with more RAM.

## Full-resolution grids: cell-centred finite volume (`--fv --capillary deck`)

The default (nodal) model stores one saturation per mesh node, shared by the cells around it. On 1-m layers
it cannot hold the sharp gas-water contact or the different saturations of neighbouring rock types, and the
full 92-layer Croft grid does not converge. With `--fv` the model stores one value per cell, like Eclipse:

    python em2ex.py -f --convert-to-m -o croft_92.e croft.grdecl
    python em2ex2moose.py --mesh croft_92.e --case CROFT_files\HM_MODEL.DATA --capillary deck --fv -o runs\croft_full_fv
    mpiexec -n 12 thermal_h2-opt -i croft_full_fv.i          (needs moose_addon_fv/, see its README)

* `--capillary deck`: the deck's Pc tables (Brooks-Corey fits, < 4 % error), scaled cell by cell like Eclipse
  (PCW/SWLPC, i.e. SWATINIT honoured), gas rel-perm end point at the cell SWL, initial water pressure =
  restart PRESSURE - PCOW: the model starts in Eclipse's own equilibrium.
* The flux is `FVPorousFlowAdvectiveFluxTPFA` (add-on): Eclipse's corner-point transmissibilities (they match
  Eclipse's TRANX/Y/Z on the full Croft grid) and gravity balanced on dipping cells. MOOSE's own FV flux
  diverges on Croft.
* Rock compressibility from ROCK with Eclipse's own formula (`ADPorousFlowPorosityRock`, add-on); `--no-rock-compressibility`
  turns it off.
* Aquifer: Eclipse's Fetkovich form - aquifer pressure at the cell depth against the cell pressure, averaged
  over each time step with the aquifer time constant (as Eclipse does), so long steps don't overdraw it.
* Croft 92 layers: the first steps converge in 2-3 Newton iterations; MUMPS needs about 5.4 GB.
* Croft 19 layers, full history (`validation/croft_fv19_vs_eclipse.png`, `..._table.txt`): 63 min on 2 cores,
  90 steps, no failed steps. Gas in place within 0.2 % of Eclipse throughout (3.376 vs 3.382 Bscf at the end);
  water influx 3.138 vs 3.091 MMstb at the end (+1.5 %), but 5-10 % high between days 60 and 400; field
  (gas) pressure within 1-2 % except around the end of production (day 365: 2440 vs 2335 psia). Whether the
  remaining mid-history difference comes from the 19-layer upscaling is best checked with the 92-layer run.
* Wells (`--well-model connection`, default with `--fv`): gas production is split over the completions by
  COMPDAT CF x gas mobility (from the last step) and each completion is limited to what it can deliver with
  BHP = 1 atm at the WELSPECS reference depth, as Eclipse's WCONHIST wells. A completion whose gas is trapped or
  gone produces nothing, so the well can no longer empty its cell (with the old fixed-rate sinks the coarse
  19-layer well cell reached Sg < 0). New outputs: `<well>_gas_rate_actual_kgps`, `<well>_cum_gas_actual_kg`,
  and `<well>_wbp`, Eclipse's WBP (CF-weighted gas pressure of the completion cells at the reference depth,
  WPAVE ... RES); `compare_with_eclipse.py` uses both. `--well-model rate` gives the old point sinks.
* Water production (`--well-water on`, default with connection wells): each completion also produces water at
  the drawdown its gas share needs, CF x rho_w krw/mu_w x dp (water mobility from the last step), as Eclipse
  produces every mobile phase; the history gas rate stays the target. Outputs `<well>_water_rate_actual_kgps`
  and `<well>_cum_water_actual_kg`; `compare_with_eclipse.py` plots them against Eclipse's WWPR/WWPT (the
  simulated water rate of the Croft report, Fig 4.6; no measured water rates exist). `--well-water off`: gas only.
* Rel-perm (`--relperm deck`, default with `--fv`): the deck's SWFN/SGFN tables with Eclipse's hysteresis
  (Croft: `SATOPTS HYSTER`, `EHYSTR 1* 3 2* KR`): gas on Killough scanning curves towards the IMBNUM tables
  (trapped gas 0.05-0.27), water on the imbibition tables scaled to IKRW = 0.55 in every cell. Needs the add-on
  `ADPorousFlowRelativePermeabilityHysteresisTable`. `--relperm corey` gives the earlier Corey fits.
  Croft 19 layers (`validation/croft_fv19_hysteresis_vs_eclipse.png`): field pressure within 3-7 psi of
  Eclipse to day 240 (Corey fits: 16-64 psi high), water influx within ~4 % (Corey: 10-40 % high).
  KNOWN LIMIT: the history wells are fixed-rate sinks. On 19 layers the single coarse completion cell is
  surrounded by water-swept cells by day ~300; with hysteresis their gas is trapped, the sink empties the
  well cell (Sg < 0) and Newton fails. (The Corey run had the same negative Sg at days 330-510 but got
  through.) A connection-based well (Eclipse CF, mobility, BHP >= 1 atm) is the fix.
* Faults (`--nnc on`, default with `--fv`): where a fault has throw the MOOSE cells on its two sides share no
  face, so without help the fault is sealed. The converter takes Eclipse's own connections for those cells
  (TRANX/Y/Z from the INIT, NNC/TRANNNC from the EGRID and INIT) and writes them to `<prefix>_nnc.csv`;
  `FVPorousFlowNNCFlux` (add-on) puts the flow on them (implicit pressure, mobility from the start of the
  step). Croft: 27,679 such pairs on 92 layers, 3,002 on 19. A fault next to CROFT_1 (I = 36|37, throw up to
  1.5 m) was sealed before: the gas west of it stayed 150-240 psi above the well side, and WBP fell 36 psi
  (day 60) to 160 psi (day 300) below Eclipse. With the connections, 19-layer hysteresis run
  (`validation/croft_fv19_nnc_wbp.png`): WBP within 1 psi of Eclipse to day 200, 3-6 psi to day 300, 21 psi
  at day 330; cell gas pressures against the Eclipse restarts (days 60-330) median within 1-5 psi (12 psi
  at day 330), was 35-95 psi. `--fault-snap` (merging nodes across small throws) is now off by default
  with `--nnc on`; `--nnc off` gives the old behaviour (snap 1 m).
* `--capillary deck` also works with the nodal model on coarse grids, but slowly; use it with `--fv`.

## Limits (same as the validated model)

* Wells are rate-controlled point sources; BHP limits are not enforced.
* Capillary pressure is zero; water held in the gas zone by capillarity in Eclipse is made immobile instead.
* Rel-perm curves are Corey fits, one set per SATNUM; Eclipse's cell-by-cell end-point scaling is not
  reproduced (initial saturations still follow the cell values).
* Rock compressibility is ignored.
* Faults with throw are no-flow (em2ex leaves those faces unconnected) in the nodal model; the `--fv` model
  connects them with Eclipse's transmissibilities (`--nnc`). Fault transmissibility multipliers (MULTFLT, MULTX
  on conforming faces) are not applied to faces MOOSE connects itself (Croft has none).

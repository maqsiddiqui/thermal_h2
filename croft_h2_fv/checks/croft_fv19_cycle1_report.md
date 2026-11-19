# Hydrogen storage input c19w.i

History input: `../hist19/fv19.i`, history output: `../hist19/fv19_exodus.e`

End-of-history state from time step 90, t = 114782400 s (1328.5 days); history end_time 114782400 s.

Aquifer 1: cumulative influx at the end of history 5.09e+08 kg (3.138 MMstb); datum pressure re-based from 19822198 Pa to 17658777 Pa (313.8 psi lower).

Mesh `/home/claude/h2fv/hist19/fv19.e` copied to `c19w.e`; 48805 cells matched to the history output (largest centroid distance 0 m); p_init and sg_init now hold the end-of-history pwater and sgas.

Initial state: 1098 gas cells (Sg > 1e-3), gas pore volume 5.847e+05 m3, pwater 16.565-22.154 MPa.

Gas in place at the end of history (history run): 7.63642e+07 kg.

Well CROFT_1: sum of open COMPDAT CF = 3.47388 cP.rb/d/psi = 9.27137e-13 m3 (deck /tmp/claude-0/-home-claude/f8745664-a354-5c1f-b624-a5c17e7f48a5/scratchpad/croft_in/CROFT_files/HM_MODEL.DATA); BHP reference depth 1984.00 m.

Injection limit: BHP <= 2865.9 psia at the reference depth (initial reservoir pressure there: 2865.9 psia; gradient 145.7 kg/m3). Withdrawal limit: BHP >= 2146.8 psia.

H2 density fit rho = p / (1.51021e+06 + 0.00831632 p), from h2_mix_pvt.csv.

Schedule: 1 cycles, injection 0.0633762 kg/s for 182.625 d, withdrawal 0.126752 kg/s for 91.3125 d; end 23668200 s.

Copied from the history input: 532 capillary-pressure objects, 532 ppss materials, 1067 other rock materials; well CROFT_1 (point 1737.731 1789.268 -1964.464); aquifer boundaries aquifer_1_Ip aquifer_1_Jm.

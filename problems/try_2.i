
#
# The mesh used below has very high resolution, so the simulation takes a long time to complete.
# Some suggested meshes of different resolution:
# nx=50, bias_x=1.2
# nx=100, bias_x=1.1
# nx=200, bias_x=1.05
# nx=400, bias_x=1.02
# nx=1000, bias_x=1.01
# nx=2000, bias_x=1.003
[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = 200
    bias_x = 1.05
    xmin = 0.1
    xmax = 5000
    ny = 1
    ymin = 0
    ymax = 11
  []
  
  [Problem]
    coord_type = RZ
  []
  
  [GlobalParams]
    PorousFlowDictator = dictator
    gravity = '0 0 0'
  []

  [UserObjects]
    [dictator]
      type = PorousFlowDictator
      porous_flow_vars = 'pgas xH2 temp'
      number_fluid_phases = 1
      number_fluid_components = 2
    []
    [pc]
      type = PorousFlowCapillaryPressureConst
      pc = 0
    []
  []
  
  [Variables]
    [pgas]
      initial_condition = 18.3e6
    []
    [xH2]
      initial_condition = 0.0
    []
    [temp]
      initial_condition = 358
    []
  []
  
  [AuxVariables]
    [rate]
    []

    #[massfrac_ph0_sp0]
    #  initial_condition = 1 # all CH4 in phase=0
    #[]
    #[massfrac_ph1_sp0]
    #  initial_condition = 0 # no CH4 in phase=1
    #[]
    #[pgas]
    #  family = MONOMIAL
    #  order = FIRST
    #[]
    #[swater]
    #  family = MONOMIAL
    #  order = FIRST
    #[]

  []
  
  [Kernels]
    [mass_gas_dot]
      type = PorousFlowMassTimeDerivative
      fluid_component = 1
      variable = pgas
      multiply_by_density = true
    []
    [flux_gas]
      type =  PorousFlowFullySaturatedAdvectiveFlux #PorousFlowAdvectiveFlux # PorousFlowFullySaturatedDarcyFlow
      fluid_component = 1
      use_displaced_mesh = false
      variable = pgas
      gravity = '0 0 0'
    []
    [mass_H2_dot]
      type = PorousFlowMassTimeDerivative
      fluid_component = 0
      variable = xH2
      multiply_by_density = true
    []
    [flux_H2]
      type = PorousFlowFullySaturatedAdvectiveFlux # PorousFlowAdvectiveFlux #PorousFlowFullySaturatedDarcyFlow
      fluid_component = 0
      use_displaced_mesh = false
      variable = xH2
      gravity = '0 0 0'
    []
    [energy_dot]
      type = PorousFlowEnergyTimeDerivative
      variable = temp
    []
    [advection]
      type = PorousFlowHeatAdvection
      use_displaced_mesh = false
      variable = temp
    []
    [conduction]
      type = PorousFlowExponentialDecay
      use_displaced_mesh = false
      variable = temp
      reference = 358
      rate = rate
    []
  []
  
  [AuxKernels]
    [rate]
      type = FunctionAux
      variable = rate
      execute_on = timestep_begin
      function = decay_rate
    []
  []
  
  [Functions]
    [decay_rate]
  # Eqn(26) of the first paper of LaForce et al.
  # Ka * (rho C)_a = 10056886.914
  # h = 11
      type = ParsedFunction
      expression = 'sqrt(10056886.914/t)/11.0'
    []
  []
  
  [FluidProperties]
    [methane]
      type = SimpleFluidProperties
      bulk_modulus = 2.27e14
      density0 = 970.0
      viscosity = 0.3394e-3
      cv = 4149.0
      cp = 4149.0
      porepressure_coefficient = 0.0
      thermal_expansion = 0
    []
    [H2]
      type = SimpleFluidProperties
      bulk_modulus = 2.27e14
      density0 = 516.48
      viscosity = 0.0393e-3
      cv = 2920.5
      cp = 2920.5
      porepressure_coefficient = 0.0
      thermal_expansion = 0
    []
    [hydrogen]
      #type = HydrogenFluidProperties # BrineFluidProperties #   HydrogenFluidProperties # BrineFluidProperties # #HydrogenFluidProperties #BrineFluidProperties
      type = SimpleFluidProperties
      bulk_modulus = 2.27e14
      density0 = 105
      viscosity = 1e-5
      cv = 2920.5
      cp = 2920.5
      porepressure_coefficient = 0.0
      thermal_expansion = 0 
    []
  []
  
  [Materials]
    [temperature]
      type = PorousFlowTemperature
      temperature = temp
    []
    [ppss]
        type = PorousFlow1PhaseFullySaturated
        porepressure = pgas
    []
    [massfrac]
      type = PorousFlowMassFraction
      mass_fraction_vars = xH2
    []
    [hydrogen]
      type =  PorousFlowSingleComponentFluid # PorousFlowMultiComponentFluid # PorousFlowSingleComponentFluid #
      fp = hydrogen
      #x = xH2
      phase = 0
    []

    [porosity_reservoir]
      type = PorousFlowPorosityConst
      porosity = 0.2
    []
    [permeability_reservoir]
      type = PorousFlowPermeabilityConst
      permeability = '2e-12 0 0  0 0 0  0 0 0'
    []
    [relperm]
     type = PorousFlowRelativePermeabilityConst
     phase = 0
    []

    [thermal_conductivity_reservoir]
      type = PorousFlowThermalConductivityIdeal
      dry_thermal_conductivity = '0 0 0  0 1.320 0  0 0 0'
      wet_thermal_conductivity = '0 0 0  0 3.083 0  0 0 0'
    []
    [internal_energy_reservoir]
      type = PorousFlowMatrixInternalEnergy
      specific_heat_capacity = 1100
      density = 2350.0
    []

  []
  
  [BCs]
    
    [cold_H2]
        type = DirichletBC
        boundary = 'left'
        variable = temp
        value = 294
    []

    [outer_temp_fixed]
      type = DirichletBC
      boundary = 'right'
      value = 358
      variable = temp
    []

    #[outer_pressure_fixed]
    #    type = DirichletBC
    #    boundary = right
    #    value = 18.3e6
    #    variable = pgas
    #[]

    [left_xH2]
        type = DirichletBC
        boundary = 'left'
        variable = xH2
        value = 1
    []

    #[left_pgas]
    #    type = DirichletBC
    #    boundary = 'left'
    #    variable = pgas
    #    value = 25e6
    #[]

    [right_pgas]
        type = DirichletBC
        boundary = 'right'
        variable = pgas
        value = 18.3e6
    []

    [H2_injection]
      type = PorousFlowSink
      boundary = left
      variable = pgas
      use_mobility = false
      use_relperm = false
      fluid_phase = 0
      flux_function = 'min(t/100.0,1)*(-2.294001475)' # 5.0E5 T/year = 15.855 kg/s, #over area of 2Pi*0.1*11
      #flux_function = -0.000125 
      #value = -0.000125
      mass_fraction_component = 1
    []

    #[right_pgas]
    #  type = PorousFlowPiecewiseLinearSink
    # boundary = 'right'
    #  variable = pgas
    #  fluid_phase = 0
    #  pt_vals = '-1E9 1E9'
    #  multipliers = '-1E9 1E9'
    # PT_shift = 18.3E6
    #  mass_fraction_component = 1
    #  use_mobility = false
    #  use_relperm = false
    #  flux_function = 1e4 #10 # 1/L
    #[]

    #[right_xH2]
    #    type = PorousFlowPiecewiseLinearSink
    #    boundary = 'right'
    #    variable = xH2
    #    fluid_phase = 0
    #   pt_vals = '-1E9 1E9'
    #    multipliers = '-1E9 1E9'
    #    PT_shift = 18.3E6
    #    mass_fraction_component = 0
    #    use_mobility = false
    #    use_relperm = false
    #    flux_function = 1e4 # 10 # 1/L
    #  []

    [outflow_a]
        type = PorousFlowOutflowBC
        boundary = right
        include_relperm = false # no need for relperm in this fully-saturated simulation
        mass_fraction_component = 0
        variable = xH2
      []
  []
  
  [Postprocessors]
    [p_bh]
      type = PointValue
      variable = pgas
      point = '0.1 0 0'
      execute_on = timestep_begin
      use_displaced_mesh = false
    []
  []
  
  [VectorPostprocessors]
    [ptsuss]
      type = LineValueSampler
      use_displaced_mesh = false
      start_point = '0.1 0 0'
      end_point = '5000 0 0'
      sort_by = x
      num_points = 50000
      outputs = csv
      variable = 'pgas temp xH2'
    []
  []
  
  [Preconditioning]
    active = 'smp'
    [smp]
      type = SMP
      full = true
      #petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
      petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
      petsc_options_value = 'gmres      asm      lu           NONZERO                   2               1E2       1E-5        500'
    []
    [mumps]
      type = SMP
      full = true
      petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
      petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -snes_rtol -snes_atol -snes_max_it'
      petsc_options_value = 'gmres      lu       mumps                         NONZERO               1E-5       1E2       50'
    []
  []
  
  [Executioner]
    type = Transient
    solve_type = NEWTON
    end_time = 1.577e+7 # 6 months #1.5768e8 #5 yrs
    #dtmax = 1e6
    [TimeStepper]
      type = IterationAdaptiveDT
      dt = 1
      growth_factor = 1.1
    []
  []
  
  [Outputs]
    print_linear_residuals = false
    sync_times = '3600 86400 2.592E6 1.577e+7' #1 hr, 1 day, 1 month, 5 yrs
    perf_graph = true
    exodus = true
    [csv]
     type = CSV
      sync_only = true
    []
  []
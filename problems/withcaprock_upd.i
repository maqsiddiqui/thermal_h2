
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
    construct_side_list_from_node_list=true
    [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 200
    bias_x = 1.025
    xmin = 0.1
    xmax = 700
    ny = 10
    ymin = 0
    ymax = 10
    
    []

  [reservoir]
    type = ParsedSubdomainMeshGenerator
    block_id = 1
    #bottom_left = '0.1 0 0'
    #top_right = '700 5 0'
    input = gen
    block_name = reservoir
    combinatorial_geometry = 'x >= 0.1 & x <= 700 & y >= 0 & y<=5'
  []

  [injection_area]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'x<0.1001'
    included_subdomains = 1
    new_sideset_name = 'injection_area'
    input = 'reservoir'
  []

  #[rename]
  #  type = RenameBlockGenerator
  #  old_block = '0 1'
  #  new_block = 'caps reservoir'
  # input = 'injection_area'
  #[]

  #[injection_nodes]
  #  input = reservoir
  #  type = ExtraNodesetGenerator
  #  new_boundary = injection_nodes
  #  coord = '0.1 0 0; 0.1 5 0'
  # []


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
      initial_condition = 14.9e6 #2170 psi
      
    []
    [xH2]
      initial_condition = 0.0
      
    []
    [temp]
      initial_condition = 331 #TVD 1470 m , 2.55C/100m
      scaling = 1E-6 # fluid enthalpy is roughly 1E6
      
    []
  []
  
  [AuxVariables]
    [rate]
    []

    #[massfrac_ph0_sp0]
    #  initial_condition = 0 # no liq H2 in phase=0
    #[]
    #[massfrac_ph1_sp0]
    #  initial_condition = 1 # all gas H2 in phase=1
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
      #gravity = '0 0 0'
      
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
      #gravity = '0 0 0'
      
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
      type = PorousFlowHeatConduction #PorousFlowExponentialDecay
      use_displaced_mesh = false
      variable = temp
      #reference = 331
      #rate = rate
      
    []
  []
  
  #[AuxKernels]
  #  [rate]
  #    type = FunctionAux
  #    variable = rate
  #    execute_on = timestep_begin
  #   function = decay_rate
  #  []
  #[]
  
  #[Functions]
  #  [decay_rate]
  # Eqn(26) of the first paper of LaForce et al.
  # Ka * (rho C)_a = 10056886.914
  # h = 11
  #    type = ParsedFunction
  #    expression = 'sqrt(10056886.914/t)/11.0'
  #  []
  #[]
  
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
    [CH4]
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
      type = HydrogenFluidProperties # BrineFluidProperties #   HydrogenFluidProperties # BrineFluidProperties # #HydrogenFluidProperties #BrineFluidProperties
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
      #fp = H2

      #x = xH2
      phase = 0
      
    []

    [porosity_reservoir]
      type = PorousFlowPorosityConst
      block = reservoir
      porosity = 0.22
    []

    [porosity_caps]
        type = PorousFlowPorosityConst
        block = 0
        porosity = 0.1
    []
    
    [permeability_reservoir]
      type = PorousFlowPermeabilityConst
      block = reservoir
      permeability = '3e-14 0 0  0 0 0  0 0 0' #3e-13 = 300 mD
    []
    
    [permeability_caps]
        type = PorousFlowPermeabilityConst
        block = 0
        permeability = '5E-16 0 0   0 5E-16 0   0 0 0'
    []  

    [relperm]
     type = PorousFlowRelativePermeabilityConst
     phase = 0
    []

    [thermal_conductivity_reservoir]
      type = PorousFlowThermalConductivityIdeal
      block = reservoir
      dry_thermal_conductivity = '1.320 0 0  0 1.320 0  0 0 1.320'
      wet_thermal_conductivity = '3.083 0 0  0 3.083 0  0 0 3.083'
    []

    [thermal_conductivity_caps]
        type = PorousFlowThermalConductivityIdeal
        block = 0
        dry_thermal_conductivity = '1.600 0 0  0 1.600 0  0 0 0'
        wet_thermal_conductivity = '4.310 0 0  0 4.310 0  0 0 0'
    []

    [internal_energy_reservoir]
      type = PorousFlowMatrixInternalEnergy
      block = reservoir
      specific_heat_capacity = 1100
      density = 2350.0
    []

    [internal_energy_caps]
        type = PorousFlowMatrixInternalEnergy
        block = 0
        specific_heat_capacity = 829
        density = 2773.4
    []

  []
  
  [BCs]
    
    [cold_H2]
        type = DirichletBC
        boundary =  'injection_area'
        variable = temp
        value = 294
    []

    [outer_temp_fixed]
      type = DirichletBC
      boundary = 'right'
      value = 331
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
        boundary = 'injection_area'
        variable = xH2
        value = 1.01
        #preset = true
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
        value = 14.9e6
    []

    #[H2_injection]
    #    type = PorousFlowSink
    #    boundary = injection_area
    #    variable = xH2
    #    #flux_function =  -0.00606 #600 ton/yr
    #    flux_function = 'min(t/100.0,1)*(-2.294001475/4)' # 5.0E5 T/year

    #[]



    [H2_injection]
      type = PorousFlowSink
      boundary = injection_area
      variable = pgas
      use_mobility = false
      use_relperm = false
      fluid_phase = 0
      #flux_function = 'min(t/100.0,1)*(-2.294001475/4)' # 5.0E5 T/year = 15.855 kg/s, #over area of 2Pi*0.1*11
      flux_function =  -0.00606 #600 ton/yr
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
        #multiply_by_density = false
      []
  []
  
  [Postprocessors]
    [p_bh]
      type = PointValue
      variable = pgas
      point = '0.1 2.5 0'
      execute_on = timestep_begin
      use_displaced_mesh = false
    []

    [injection_area]
      type = AreaPostprocessor
      boundary = injection_area
      execute_on = initial
    []
    [injected_H2]
      type = PorousFlowFluidMass
      fluid_component = 0
    []

  []
  
  [VectorPostprocessors]
    [ptsuss]
      type = LineValueSampler
      use_displaced_mesh = false
      start_point = '0.1 2.5 0'
      end_point = '700 2.5 0'
      sort_by = x
      num_points = 7000
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
      petsc_options_value = 'gmres      asm      lu           NONZERO                   2               1E-5       1E-5        500'
    []
    [mumps]
      type = SMP
      full = true
      petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
      petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -snes_rtol -snes_atol -snes_max_it'
      petsc_options_value = 'gmres      lu       mumps                         NONZERO               1E-5      1E-5       50'
    []
  []
  
  [Executioner]
    type = Transient
    solve_type = NEWTON
    end_time = 1.577e+7 # 6 months #1.5768e8 #5 yrs
    dtmax = 1e6
    [TimeStepper]
      type = IterationAdaptiveDT
      dt = 1000
      growth_factor = 2
    []
    nl_abs_tol = 1e-8
  []

  [Debug]
    show_var_residual_norms = true
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
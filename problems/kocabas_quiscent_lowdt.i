########
#### THIS CODE SOLVES THE ADVONIN 1965. 2D RADIAL THERMAL FLOW PROBLEM :
########

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
    bias_x = 1.05
    xmin = 0.1
    xmax = 1000
    ny = 150
    ymin = 0
    ymax = 300
    #rz_coord_axis = Y
    #coord_type = RZ
    []

    [reservoir]
        type = ParsedSubdomainMeshGenerator
        block_id = 1
        #bottom_left = '0.1 0 0'
        #top_right = '700 5 0'
        input = gen
        block_name = reservoir
        combinatorial_geometry = 'x >= 0.1 & x <= 1000 & y >= 0 & y<=50'
    []
      
    [injection_area]
        type = ParsedGenerateSideset
        combinatorial_geometry = 'x<0.1001'
        included_subdomains = 1
        new_sideset_name = 'injection_area'
        input = 'reservoir'
      []

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
      porous_flow_vars = 'pgas temp'
      number_fluid_phases = 1
      number_fluid_components = 1
    []
    [pc]
      type = PorousFlowCapillaryPressureConst
      pc = 0
    []
  []
  
  [Variables]
    [pgas]
      initial_condition = 5e6 #2170 psi
    []
    [temp]
      initial_condition = 443.15 #TVD 1470 m , 2.55C/100m
      scaling = 1E-6 # fluid enthalpy is roughly 1E6
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
      fluid_component = 0
      variable = pgas
      multiply_by_density = true
      
    []
    [flux_gas]
      type =  PorousFlowAdvectiveFlux # PorousFlowFullySaturatedDarcyFlow
      fluid_component = 0
      use_displaced_mesh = false
      variable = pgas
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
  

  [Functions]
    [injection_rate]
      type = ParsedFunction
      #symbol_values = injection_area
      #symbol_names = area
      #expression = '-0.5*9.052/area' #0.01 m3/s * rhow
      value = -0.5*9052400000/(2*3.14*0.1*100*1e9)

    []
  []

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
    [water]
      type = SimpleFluidProperties
      bulk_modulus = 6.7e10 #2.27e14
      density0 = 905.24
      viscosity = 0.1658e-3
      cv = 3456.2
      cp = 4336.88
      porepressure_coefficient = 0.0
      thermal_expansion = 0
      thermal_conductivity = 0.6817
    []
    #[hydrogen]
    #  type = HydrogenFluidProperties # BrineFluidProperties #   HydrogenFluidProperties # BrineFluidProperties # #HydrogenFluidProperties #BrineFluidProperties
    #[]
  []
  
  [Materials]
    [temperature]
      type = PorousFlowTemperature
      temperature = temp
      #at_nodes = true
    []
    [ppss]
        type = PorousFlow1PhaseFullySaturated
        porepressure = pgas
    []

    [massfrac]
        type = PorousFlowMassFraction
        mass_fraction_vars = ''
    []

    [water]
      type =  PorousFlowSingleComponentFluid # PorousFlowMultiComponentFluid # PorousFlowSingleComponentFluid #
      fp = water
      #fp = H2
      #x = xH2
      phase = 0
      #at_nodes = true
    []

    [porosity_reservoir]
      type = PorousFlowPorosityConst
      porosity = 0.20
      #at_nodes = true
      block = reservoir
    []

    [porosity_caps]
        type = PorousFlowPorosityConst
        block = 0
        porosity = 0.1
    []

    [permeability_reservoir]
      type = PorousFlowPermeabilityConst
      #permeability = '1.8e-11 0 0  0 1.8e-11 0  0 0 1.8e-11' #3e-13 = 300 mD
      permeability = '6.6646e-13 0 0  0 6.6646e-13 0  0 0 6.6646e-13'
      block = reservoir
    []

    [permeability_caps]
        type = PorousFlowPermeabilityConst
        block = 0
        permeability = '5E-18 0 0   0 5E-18 0   0 0 0'
    [] 

    [relperm]
     type = PorousFlowRelativePermeabilityConst
     phase = 0
     #at_nodes = true
    []

    [thermal_conductivity_reservoir]
      type = PorousFlowThermalConductivityIdeal
     
      dry_thermal_conductivity = '0 0 0  0 20 0  0 0 0'
      #wet_thermal_conductivity = '3.083 0 0  0 0 0  0 0 3.083'
      #wet_thermal_conductivity = '0 0 0  0 0 0  0 0 3.083'
      block = reservoir
    []

    [thermal_conductivity_caps]
        type = PorousFlowThermalConductivityIdeal
        block = 0
        dry_thermal_conductivity = '0 0 0  0 20 0  0 0 0'
        #wet_thermal_conductivity = '0 0 0  0 4.310 0  0 0 0'
    []

    [internal_energy_reservoir]
      type = PorousFlowMatrixInternalEnergy
      specific_heat_capacity = 739.95
      density = 2898.68
      block = reservoir
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
        boundary =  injection_area
        variable = temp
        value = 433.15
    []

    [outer_temp_fixed]
      type = DirichletBC
      boundary = 'top right'
      value = 443.15
      variable = temp
    []

    #[outer_pressure_fixed]
    #    type = DirichletBC
    #    boundary = right
    #    value = 18.3e6
    #    variable = pgas
    #[]

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
        value = 5e6
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
      #fluid_phase = 0
      #flux_function = 'min(t/100.0,1)*(-2.294001475/4)' # 5.0E5 T/year = 15.855 kg/s, #over area of 2Pi*0.1*11
      #flux_function =  -0.00606 #600 ton/yr
      flux_function =  injection_rate
      #value = -0.000125
      #mass_fraction_component = 1
    []

    [H2_injection_stop]
        type = PorousFlowSink
        boundary = injection_area
        variable = pgas
        use_mobility = false
        use_relperm = false
        #fluid_phase = 0
        #flux_function = 'min(t/100.0,1)*(-2.294001475/4)' # 5.0E5 T/year = 15.855 kg/s, #over area of 2Pi*0.1*11
        #flux_function =  -0.00606 #600 ton/yr
        flux_function =  0
        #value = -0.000125
        #mass_fraction_component = 1
      []

    #[outflow_a]
    #    type = PorousFlowOutflowBC
    #    boundary = right
    #    include_relperm = false # no need for relperm in this fully-saturated simulation
    #    #mass_fraction_component = 0
    #   variable = pgas
    #  []
[]

[Controls]
    [inject]
      type = TimePeriod
      enable_objects = 'BCs::H2_injection'
      disable_objects = 'BCs::H2_injection_stop'
      start_time = 0
      end_time = 1e9
      set_sync_times = true
      execute_on = 'initial timestep_begin'
      implicit = false
    []

    [inject_stop]
        type = TimePeriod
        enable_objects = 'BCs::H2_injection_stop'
        disable_objects = 'BCs::H2_injection'
        start_time = 1e9
        end_time = 2e9
        set_sync_times = true
        execute_on = 'initial timestep_begin'
        implicit = false
    []

    [temp]
        type = TimePeriod
        enable_objects = 'BCs::cold_H2'
        start_time = 0
        end_time = 1e9
        set_sync_times = true
        execute_on = 'initial timestep_begin'
        implicit = false
      []
[]
  
[Postprocessors]

    [injection_area]
        type = AreaPostprocessor
        boundary = injection_area
        execute_on = initial
    []

    #[injection_area]
    #  type = AreaPostprocessor
    #  boundary = injection_area
    #  execute_on = initial
    #[]
    #[injected_H2]
    #  type = PorousFlowFluidMass
    #  fluid_component = 0
    #[]

[]
  
  [VectorPostprocessors]
    [ptsussq1HC] #High caprock
      type = LineValueSampler
      use_displaced_mesh = false
      start_point = '0.1 25 0'
      end_point = '1000 25 0'
      sort_by = x
      num_points = 1000
      outputs = csv
      variable = 'pgas temp'
    []

    [ptsussq2HC]
            type = LineValueSampler
            use_displaced_mesh = false
            start_point = '0.1 50 0'
            end_point = '1000 50 0'
            sort_by = x
            num_points = 1000
            outputs = csv
            variable = 'pgas temp'
    []

    [ptsussq3HC]
            type = LineValueSampler
            use_displaced_mesh = false
            start_point = '0.1 75 0'
            end_point = '1000 75 0'
            sort_by = x
            num_points = 1000
            outputs = csv
            variable = 'pgas temp'
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
    end_time = 2e9
    dtmax = 1e5
    [TimeStepper]
      type = IterationAdaptiveDT
      dt = 100
      growth_factor = 2
    []
    nl_abs_tol = 1e-8
  []

  [Debug]
    show_var_residual_norms = true
  []
  
  [Outputs]
    print_linear_residuals = false
    sync_times = '1.1e9 2e9' #after injection 
    perf_graph = true
    exodus = true
    [csv]
     type = CSV
      sync_only = true
    []
  []
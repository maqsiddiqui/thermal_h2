# inject_fluid_mass = 180e3 # kg

inject_time = 1.577e7 # # 6 months in seconds
store_time =  0 # seconds
produce_time = 1.577e7 #1.577e7 # # 6 months in seconds
rest_time = 0 # days
num_cycles = 1 # Currently needs to be <= 10

cycle_length = '${fparse inject_time + store_time + produce_time + rest_time}'

end_simulation = '${fparse cycle_length * num_cycles}'

start_injection1 = 0

end_injection1 = '${fparse start_injection1 + inject_time}'

start_production1 = '${fparse end_injection1 + store_time}'

end_production1 = '${fparse start_production1 + produce_time}'


synctimes = '${start_injection1} ${end_injection1} ${end_production1}'


[Mesh]
    construct_side_list_from_node_list=true
    [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 500
    # bias_x = 1.05
    xmin = 0.1
    xmax = 1000
    ny = 1
    ymin = 0
    ymax = 50
    # rz_coord_axis = Y
    # coord_type = RZ
    []

    [injection_area]
        type = ParsedGenerateSideset
        combinatorial_geometry = 'x<0.1001'
        # included_subdomains = 1
        new_sideset_name = 'injection_area'
        input = 'gen'
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
      porous_flow_vars = 'porepressure cgH2 cgCH4'
      number_fluid_phases = 1
      number_fluid_components = 3
    []
    [pc]
      type = PorousFlowCapillaryPressureConst
      pc = 0
    []

    [fluid_advective_flux]
        type = PorousFlowAdvectiveFluxCalculatorSaturated
        flux_limiter_type = superbee
    []

[]
  
[Variables]
    [porepressure]

    []

    [cgH2]
      initial_condition = 0 #2170 psi
      scaling = 1e4
    []
    [cgCH4]
      initial_condition = 0 #
      scaling = 1e4
    []
[]

[ICs]

    [porepressure]
        type = FunctionIC
        variable = porepressure
        function = 14.5e6-x
    []
[]
  
[AuxVariables]
    [rate]
    []
[]

  
[Kernels]

    [mass_gas_dot]
        type = PorousFlowMassTimeDerivative
        fluid_component = 0
        variable = porepressure
        # multiply_by_density = true
    []

    [flux_gas]
        type = PorousFlowFullySaturatedDarcyFlow
        fluid_component = 0
        variable = porepressure
    []

    [diff_gas]
        type = PorousFlowDispersiveFlux
        fluid_component = 0
        variable = porepressure
        disp_trans = 0
        disp_long = 0
    []

    [cgH2_dt]
        type = PorousFlowMassTimeDerivative
        fluid_component = 1
        variable = cgH2
        #multiply_by_density = true
    []
    [cgH2_diff]
        type = PorousFlowDispersiveFlux
        fluid_component = 1
        variable = cgH2
        disp_long = 0
        disp_trans = 0
    []
    [cgH2_conv]
        type =  PorousFlowFullySaturatedDarcyFlow
        fluid_component = 1
        variable = cgH2
    []

    [cgH2_reaction]
        type = CoefReaction
        coefficient = -9.84e-8
        variable = cgH2
    []

    [cgCH4_dt]
        type = PorousFlowMassTimeDerivative
        fluid_component = 2
        variable = cgCH4
       # multiply_by_density = true
    []
    [cgCH4_diff]
        type = PorousFlowDispersiveFlux
        fluid_component = 2
        variable = cgCH4
        disp_long = 0
        disp_trans = 0

    []
    [cgCH4_conv]
        type = PorousFlowFullySaturatedDarcyFlow
        fluid_component = 2
        variable = cgCH4
    []

    [cgCH4_reaction]
        type = CoupledForce
        coef = -2.46e-8 #
        v = cgH2
        variable = cgCH4
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
      symbol_values = injection_area
      symbol_names = area
      expression = '-0.5*9.052/area' #0.01 m3/s * rhow
    []

    [production_rate]
        type = ParsedFunction
        symbol_values = injection_area
        symbol_names = area
        expression = '0.5*9.052/area' #0.01 m3/s * rhow
    []

    [inject]
        type = ParsedFunction
        expression = 'if(t >= ${start_injection1} & t < ${end_injection1}, 1,0)'
    
    []

    [produce]
        type = ParsedFunction
        expression = 'if(t >= ${start_production1} & t < ${end_production1}, 1,0)'
    
    []

[]

[Controls]
    [inject_h2_bc]
      type = ConditionalFunctionEnableControl
      enable_objects = 'BCs::cgH2_in'
      conditional_function = inject
      implicit = false
      execute_on = 'initial timestep_begin'
    []
    [inject_h2]
        type = ConditionalFunctionEnableControl
        enable_objects = 'BCs::H2_injection'
        conditional_function = inject
        implicit = false
        execute_on = 'initial timestep_begin'
    []
    [produce_h2]
        type = ConditionalFunctionEnableControl
        enable_objects = 'BCs::H2_production BCs::CH4_production BCs::cgH2_production BCs::outflow_pp BCs::outflow_h2 BCs::outflow_ch4'
        conditional_function = produce
        implicit = false
        execute_on = 'initial timestep_begin'
    []

    [inject_ch4_bc]
      type = ConditionalFunctionEnableControl
      enable_objects = 'BCs::cgCH4_in'
      conditional_function = inject
      implicit = false
      execute_on = 'initial timestep_begin'
    []

[]
  
[FluidProperties]
    [water]
        type=SimpleFluidProperties
        bulk_modulus = 14.5e6
        density0 = 9.55 # 905.24
        viscosity = 9.93e-6
        cv = 14657.2
        cp = 14657.2
        porepressure_coefficient = 0.0
        thermal_expansion = 0
        thermal_conductivity = 0.22
    []
    #[hydrogen]
    #  type = HydrogenFluidProperties # BrineFluidProperties #   HydrogenFluidProperties # BrineFluidProperties # #HydrogenFluidProperties #BrineFluidProperties
    #[]
[]
  
[Materials]
    [temperature]
      type = PorousFlowTemperature
      temperature = 350
      #at_nodes = true
    []
    [ppss]
        type = PorousFlow1PhaseFullySaturated
        porepressure = porepressure
    []

    [massfrac]
        type = PorousFlowMassFraction
        mass_fraction_vars = 'cgH2 cgCH4'
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
      porosity = 0.18
      #at_nodes = true
    []

    [permeability_reservoir]
      type = PorousFlowPermeabilityConst
      #permeability = '1.8e-11 0 0  0 1.8e-11 0  0 0 1.8e-11' #3e-13 = 300 mD
      permeability = '6.6646e-13 0 0  0 6.6646e-13 0  0 0 6.6646e-13'
    []

    [relperm]
     type = PorousFlowRelativePermeabilityConst
     phase = 0
     #at_nodes = true
    []

    [diff]
        type = PorousFlowDiffusivityConst
        diffusion_coeff = '4.4e-5 4.4e-5 4.4e-5'
        tortuosity = 1
    []

  []
  
[BCs]

    [cgH2_in]
        type = DirichletBC
        variable = cgH2
        boundary = injection_area
        value = 1
    []
  
    [cgCH4_in]
          type = DirichletBC
          variable = cgCH4
          boundary = injection_area
          value = 0
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
      variable = porepressure
      use_mobility = false
      use_relperm = false
      #fluid_phase = 0
      #flux_function = 'min(t/100.0,1)*(-2.294001475/4)' # 5.0E5 T/year = 15.855 kg/s, #over area of 2Pi*0.1*11
      #flux_function =  -0.00606 #600 ton/yr
      flux_function =  injection_rate
      #value = -0.000125
      #mass_fraction_component = 1
    []

    [H2_production]
        type = PorousFlowSink
        boundary = injection_area
        variable = porepressure
        use_mobility = false
        use_relperm = false
        fluid_phase = 0
        #flux_function = 'min(t/100.0,1)*(-2.294001475/4)' # 5.0E5 T/year = 15.855 kg/s, #over area of 2Pi*0.1*11
        #flux_function =  -0.00606 #600 ton/yr
        flux_function =  production_rate
        #value = -0.000125
        mass_fraction_component = 0
    []

    [CH4_production]
        type = PorousFlowSink
        boundary = injection_area
        variable = cgCH4
        use_mobility = false
        use_relperm = false
        flux_function =  production_rate
        fluid_phase = 0
        mass_fraction_component = 2
    []

    [cgH2_production]
        type = PorousFlowSink
        boundary = injection_area
        variable = cgH2
        use_mobility = false
        use_relperm = false
        fluid_phase = 0
        flux_function =  production_rate
        mass_fraction_component = 1
    []

    # [inject_right]
    #     type = PorousFlowOutflowBC
    #     variable = porepressure
    #     boundary = right
    #     mass_fraction_component = 0
    #     # flux_function = injection_rate_value
    # []

    [outflow_pp]
       type = PorousFlowOutflowBC
       boundary = left
       include_relperm = false # no need for relperm in this fully-saturated simulation
       mass_fraction_component = 0
       variable = porepressure
     []

    [outflow_h2]
        type = PorousFlowOutflowBC
        boundary = left
        include_relperm = false # no need for relperm in this fully-saturated simulation
        mass_fraction_component = 1
        variable = cgH2
        
    []

    [outflow_ch4]
        type = PorousFlowOutflowBC
        boundary = left
        include_relperm = false # no need for relperm in this fully-saturated simulation
        mass_fraction_component = 2
        variable = cgCH4
        
    []

[]
  
[Postprocessors]

    [injection_area]
        type = AreaPostprocessor
        boundary = injection_area
        execute_on = initial
    []

[]
  
[VectorPostprocessors]
    [ptsuss]
      type = LineValueSampler
      use_displaced_mesh = false
      start_point = '0.1 25 0'
      end_point = '1000 25 0'
      sort_by = x
      num_points = 1000
      outputs = csv
      variable = 'cgH2 cgCH4'
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
    end_time = ${end_simulation}
    dtmax = 1e4
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
    sync_times = ${synctimes} #1 hr, 1 day, 1 month, 5 yrs
    perf_graph = true
    exodus = true
    [csv]
     type = CSV
      sync_only = true
    []
[]
inject_fluid_mass = 180e3 # kg
produce_fluid_mass = ${inject_fluid_mass} # kg
#true_screen_area = 1 # 100 #m2

inject_time = 1.577e7 # # 6 months in seconds
store_time =  0 # seconds
produce_time = 1.577e7 # # 6 months in seconds
rest_time = 0 # days
num_cycles = 5 # Currently needs to be <= 10

cycle_length = '${fparse inject_time + store_time + produce_time + rest_time}'

end_simulation = '${fparse cycle_length * num_cycles}'

start_injection1 = 0
start_injection2 = ${cycle_length}
start_injection3 = '${fparse cycle_length * 2}'
start_injection4 = '${fparse cycle_length * 3}'
start_injection5 = '${fparse cycle_length * 4}'

end_injection1 = '${fparse start_injection1 + inject_time}'
end_injection2 = '${fparse start_injection2 + inject_time}'
end_injection3 = '${fparse start_injection3 + inject_time}'
end_injection4 = '${fparse start_injection4 + inject_time}'
end_injection5 = '${fparse start_injection5 + inject_time}'


start_production1 = '${fparse end_injection1 + store_time}'
start_production2 = '${fparse end_injection2 + store_time}'
start_production3 = '${fparse end_injection3 + store_time}'
start_production4 = '${fparse end_injection4 + store_time}'
start_production5 = '${fparse end_injection5 + store_time}'

end_production1 = '${fparse start_production1 + produce_time}'
end_production2 = '${fparse start_production2 + produce_time}'
end_production3 = '${fparse start_production3 + produce_time}'
end_production4 = '${fparse start_production4 + produce_time}'
end_production5 = '${fparse start_production5 + produce_time}'


# plot_time_inj1 = 7.884e6 # 3 months in secs
# plot_time_prod1 = 2.365e7 # 9 months in secs
 synctimes = '${start_injection1} ${end_injection1} ${end_production1} 
              ${end_injection2} ${end_production2}
              ${end_injection3} ${end_production3}
              ${end_injection4} ${end_production4}
              ${end_injection5} ${end_production5}'


#end_simulation = 1.577e7
#synctimes = '2.628e6 7.884e6 1.577e7' # 1 mo, 3 mo, 6 mo

[Mesh]
    type = GeneratedMesh
    dim = 1
    xmin = 0
    xmax = 1000
    nx = 500
[]

[GlobalParams]
    PorousFlowDictator = dictator
    gravity = '0 0 0'
[]

[Variables]
    [porepressure]
        # order = FIRST
        # family = LAGRANGE
    []
    [cgH2]
    #   order = FIRST
    #   family = LAGRANGE
      scaling = 1e4
    []
    [cgCH4]
        # order = FIRST
        # family = LAGRANGE
        scaling = 1e4
    []
[]

# [AuxVariables]
#     [porepressure]
#       order = FIRST
#       family = LAGRANGE
#     []
# []


[ICs]
    [cgH2]
      type = ConstantIC
      variable = cgH2
      value = 0
    []
    [cgCH4]
      type = ConstantIC
      variable = cgCH4
      value = 0
    []
    [porepressure]
        type = FunctionIC
        variable = porepressure
        function = 14.5e6-x
    []
[]


[Kernels]
    [mass_gas_dot]
        type = PorousFlowMassTimeDerivative
        fluid_component = 1
        variable = porepressure
        # multiply_by_density = true
    []

    [flux_gas]
        type = PorousFlowFullySaturatedDarcyFlow
        fluid_component = 1
        variable = porepressure
    []
    [diff_H2]
        type = PorousFlowDispersiveFlux
        fluid_component = 1
        variable = porepressure
        disp_trans = 0
        disp_long = 0
    []

    [cgH2_dt]
        type = PorousFlowMassTimeDerivative
        fluid_component = 0
        variable = cgH2
        multiply_by_density = true
    []
    [cgH2_diff]
        type = PorousFlowDispersiveFlux
        fluid_component = 0
        variable = cgH2
        disp_long = 0
        disp_trans = 0
    []
    [cgH2_conv]
        type =  PorousFlowFullySaturatedDarcyFlow
        fluid_component = 0
        variable = cgH2

    []

    [cgH2_reaction]
        type = CoefReaction
        coefficient = 9.84e-8
        variable = cgH2
    []

    [cgCH4_dt]
        type = PorousFlowMassTimeDerivative
        fluid_component = 2
        variable = cgCH4
        multiply_by_density = true
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
        coef = 2.46e-8 #
        v = cgH2
        variable = cgCH4
    []

[]

[BCs]

    # [outer_pressure_fixed]
    #     type = DirichletBC
    #     boundary = left
    #     value = 14.5e6
    #     variable = porepressure
    # []

    # [constant_outer_porepressure]
    #     type = DirichletBC
    #     variable = porepressure
    #     value = 0
    #     boundary = right
    # []

    [cgH2_in]
      type = DirichletBC
      variable = cgH2
      boundary = left
      value = 1
    []

    [cgCH4_in]
        type = DirichletBC
        variable = cgCH4
        boundary = left
        value = 0
    []

    # [cgH2_out]
    #     type = PorousFlowSink
    #     variable = cgH2
    #     boundary = left
    #     flux_function = production_rate_value
    #     fluid_phase = 0
    #     use_enthalpy = true
    #     #save_in = heat_flux_out
    # []
    # [cgCH4_out]
    #     type = PorousFlowSink
    #     variable = cgCH4
    #     boundary = left
    #     flux_function = production_rate_value
    #     fluid_phase = 0
    #     use_enthalpy = true
    #     #save_in = heat_flux_out
    # []

    [inject_fluid]
        type = PorousFlowSink
        variable = porepressure
        boundary = left
        flux_function = injection_rate_value
    []

    [produce_fluid]
        type = PorousFlowSink
        variable = porepressure
        boundary = left
        flux_function = production_rate_value
    []

    [inject_right]
        type = PorousFlowOutflowBC
        variable = porepressure
        boundary = right
        mass_fraction_component = 0
        # flux_function = injection_rate_value
    []


    # [inject_fluidm]
    #     type = PorousFlowSink
    #     variable = cgCH4
    #     boundary = left
    #     flux_function = injection_rate_value
    # []

    # [produce_fluidm]
    #     type = PorousFlowSink
    #     variable = cgCH4
    #     boundary = left
    #     flux_function = production_rate_value
    # []

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
        enable_objects = 'BCs::inject_fluid'
        conditional_function = inject
        implicit = false
        execute_on = 'initial timestep_begin'
    []
    [produce_h2]
        type = ConditionalFunctionEnableControl
        enable_objects = 'BCs::produce_fluid'
        conditional_function = produce
        implicit = false
        execute_on = 'initial timestep_begin'
    []
    # [prod_h2_kern]
    #     type = ConditionalFunctionEnableControl
    #     enable_objects = 'Kernels::cgH2_ex'
    #     conditional_function = produce
    #     implicit = false
    #     execute_on = 'initial timestep_begin'
    # []

    [inject_ch4_bc]
      type = ConditionalFunctionEnableControl
      enable_objects = 'BCs::cgCH4_in'
      conditional_function = inject
      implicit = false
      execute_on = 'initial timestep_begin'
    []
    # [prod_ch4_kern]
    #     type = ConditionalFunctionEnableControl
    #     enable_objects = 'Kernels::cgCH4_ex'
    #     conditional_function = produce
    #     implicit = false
    #     execute_on = 'initial timestep_begin'
    # []

[]

[Functions]

  [inject]
    type = ParsedFunction
    expression = 'if(t >= ${start_injection1} & t < ${end_injection1}, 1,
                  if(t >= ${start_injection2} & t < ${end_injection2}, 1,
                  if(t >= ${start_injection3} & t < ${end_injection3}, 1,
                  if(t >= ${start_injection4} & t < ${end_injection4}, 1,
                  if(t >= ${start_injection5} & t < ${end_injection5}, 1,0)))))'

  []
  [produce]
    type = ParsedFunction
    expression = 'if(t >= ${start_production1} & t < ${end_production1}, 1,
                  if(t >= ${start_production2} & t < ${end_production2}, 1,
                  if(t >= ${start_production3} & t < ${end_production3}, 1,
                  if(t >= ${start_production4} & t < ${end_production4}, 1,
                  if(t >= ${start_production5} & t < ${end_production5}, 1,0)))))'

  []

  [injection_rate_value]
    type = ParsedFunction
    symbol_names = true_screen_area
    symbol_values = injection_area
    expression = '-${inject_fluid_mass}/(true_screen_area * ${inject_time})'
  []
  [production_rate_value]
    type = ParsedFunction
    symbol_names = true_screen_area
    symbol_values = injection_area
    expression = '${produce_fluid_mass}/(true_screen_area * ${produce_time})'
  []

[]

[Postprocessors]
  [injection_area]
    type = AreaPostprocessor
    boundary = left
    execute_on = initial
  []
[]


[FluidProperties]
    [true_water]
      type = SimpleFluidProperties  # 
      bulk_modulus = 14.5e6
      density0 = 9.55 # 905.24
      viscosity = 9.93e-6
      cv = 14657.2
      cp = 14657.2
      porepressure_coefficient = 0.0
      thermal_expansion = 0
      thermal_conductivity = 0.22
    []

[]
[Materials]
    # [porous_inj]
    #     type = GenericConstantMaterial
    #     prop_names = 'diffusivity conductivity porosity'
    #     prop_values = '4.4e-5 3.27e-5 1'
    # []
    [temperature]
        type = PorousFlowTemperature
        temperature = 68
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
      fp = true_water
      #fp = H2
      #x = xH2
      phase = 0
      temperature_unit = Celsius
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
      permeability = '6.6646e-13 0 0  0 0 0  0 0 0'
    []

    [relperm]
     type = PorousFlowRelativePermeabilityConst
     phase = 0
     #at_nodes = true
    []

    [diff]
        type = PorousFlowDiffusivityConst
        diffusion_coeff = '4.4e-5 4.4e-5 4.4e-5'
        tortuosity =1
    []
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
  
[VectorPostprocessors]
    [lineplot]
      type = LineValueSampler
      use_displaced_mesh = false
      start_point = '0 0 0'
      end_point = '1000 0 0'
      sort_by = x
      num_points = 200
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
      sync_times = ${synctimes}
      perf_graph = true
      exodus = true
    [csv]
      type = CSV
      sync_only = true
    []
[]
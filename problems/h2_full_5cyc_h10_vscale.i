## This code solves temperature changes during H2 injection-withdrawal with no intermediary storage

inject_fluid_mass = 3.6e6 # kg
produce_fluid_mass = ${inject_fluid_mass} # kg
#inject_temp = 318.15 # degC

inject_time2 = 1.577e7
inject_time = '${fparse 5*inject_time2}' # secs
store_time =  0 # 
produce_time2 = 1.577e7 # # days
produce_time = '${fparse 5*produce_time2}' # secs
rest_time = 0 # days
num_cycles = 5 # Currently needs to be <= 10

cycle_length = '${fparse inject_time + store_time + produce_time + rest_time}'

end_simulation = '${fparse cycle_length * num_cycles}'

# Note: I have setup 10 cycles but you can set num_cycles less than 10.
start_injection1 = 0
start_injection2 = ${cycle_length}
start_injection3 = '${fparse cycle_length * 2}'
start_injection4 = '${fparse cycle_length * 3}'
start_injection5 = '${fparse cycle_length * 4}'
# start_injection6 = '${fparse cycle_length * 5}'
# start_injection7 = '${fparse cycle_length * 6}'
# start_injection8 = '${fparse cycle_length * 7}'
# start_injection9 = '${fparse cycle_length * 8}'
# start_injection10 = '${fparse cycle_length * 9}'

end_injection1 = '${fparse start_injection1 + inject_time}'
end_injection2 = '${fparse start_injection2 + inject_time}'
end_injection3 = '${fparse start_injection3 + inject_time}'
end_injection4 = '${fparse start_injection4 + inject_time}'
end_injection5 = '${fparse start_injection5 + inject_time}'
# end_injection6 = '${fparse start_injection6 + inject_time}'
# end_injection7 = '${fparse start_injection7 + inject_time}'
# end_injection8 = '${fparse start_injection8 + inject_time}'
# end_injection9 = '${fparse start_injection9 + inject_time}'
# end_injection10 = '${fparse start_injection10 + inject_time}'

start_production1 = '${fparse end_injection1 + store_time}'
start_production2 = '${fparse end_injection2 + store_time}'
start_production3 = '${fparse end_injection3 + store_time}'
start_production4 = '${fparse end_injection4 + store_time}'
start_production5 = '${fparse end_injection5 + store_time}'
# start_production6 = '${fparse end_injection6 + store_time}'
# start_production7 = '${fparse end_injection7 + store_time}'
# start_production8 = '${fparse end_injection8 + store_time}'
# start_production9 = '${fparse end_injection9 + store_time}'
# start_production10 = '${fparse end_injection10 + store_time}'

end_production1 = '${fparse start_production1 + produce_time}'
end_production2 = '${fparse start_production2 + produce_time}'
end_production3 = '${fparse start_production3 + produce_time}'
end_production4 = '${fparse start_production4 + produce_time}'
end_production5 = '${fparse start_production5 + produce_time}'
# end_production6 = '${fparse start_production6 + produce_time}'
# end_production7 = '${fparse start_production7 + produce_time}'
# end_production8 = '${fparse start_production8 + produce_time}'
# end_production9 = '${fparse start_production9 + produce_time}'
# end_production10 = '${fparse start_production10 + produce_time}'

# plot_time_inj1 = 7.884e6 # 3 months in secs
# plot_time_prod1 = 2.365e7 # 9 months in secs
 synctimes = '${start_injection1} ${end_injection1} ${end_production1} 
              ${end_injection2} ${end_production2}
              ${end_injection3} ${end_production3}
              ${end_injection4} ${end_production4}
              ${end_injection5} ${end_production5}'
            #  ${start_injection6} ${end_injection6} ${start_production6} ${end_production6}
            #  ${start_injection7} ${end_injection7} ${start_production7} ${end_production7}
            #  ${start_injection8} ${end_injection8} ${start_production8} ${end_production8}
            #  ${start_injection9} ${end_injection9} ${start_production9} ${end_production9}
            #  ${start_injection10} ${end_injection10} ${start_production10} ${end_production10}'


aq_thickness = 10
# injection region top and bottom (m).  Note, the mesh is created with the aquifer in y = (-0.5 * aq_thickness, 0.5 * aq_thickness), irrespective of depth (depth only sets the insitu porepressure and temperature)
screen_top = '${fparse 0.5 * aq_thickness}'
screen_bottom = '${fparse -0.5 * aq_thickness}'

#####################################
half_aq_thickness = '${fparse aq_thickness * 0.5}'

[Mesh]
    construct_side_list_from_node_list=true
    [mesh]
        type = CartesianMeshGenerator
        dim = 2
        # xmin = 0.1
        dx = '1 2 12 15 10 10 50'
        ix = '10 10 50 30 10 10 10'
        dy = '25 30 25'
        iy = '10 30 10'
    []

    [shift]
        type = TransformGenerator
        input = mesh
        transform = translate
        vector_value = '0.1 -40 0'
    []

    [subdomains]
        type = ParsedSubdomainMeshGenerator
        input = shift
        combinatorial_geometry = 'y >= -${half_aq_thickness} & y <= ${half_aq_thickness}'
        block_id = 1
    []

    [injection_area]
        type = ParsedGenerateSideset #ParsedSubdomainMeshGenerator
        combinatorial_geometry = 'x<=0.1001 & y >= ${screen_bottom} & y <= ${screen_top}'
        included_subdomains = 1
        new_sideset_name = 'injection_area'
        # block_id = 1
        input = 'subdomains'
    []

    [rename]
        type = RenameBlockGenerator
        old_block = '0 1'
        new_block = 'caps aquifer'
        input = injection_area
    []



      coord_type = RZ
    
[]
  

[GlobalParams]
    PorousFlowDictator = dictator
    gravity = '0 -9.81 0'
[]

[Variables]
    [porepressure]
      initial_condition = 5e6 #2170 psi
    []
    [temperature]
      initial_condition = 68 #341.15 #68 #341.15 #TVD 1470 m , 2.55C/100m
      scaling = 1E-6 # fluid enthalpy is roughly 1E6
    []
[]

# [PorousFlowFullySaturated]
#     coupling_type = ThermoHydro
#     porepressure = porepressure
#     temperature = temperature
#     fp = true_water
#     stabilization = KT
#     flux_limiter_type = ${flux_limiter}
#     use_displaced_mesh = false
#     temperature_unit = Celsius
#     pressure_unit = Pa
#     time_unit = days
# []

[BCs]
    [outer_boundary_porepressure]
      type = DirichletBC
      preset = true
      variable = porepressure
    #   function = insitu_pressure
      value = 5e6
      boundary = 'bottom right top'
    []
    [outer_boundary_temperature]
      type = DirichletBC
      preset = true
      variable = temperature
    #   function = insitu_temperature
      value = 68 #341.15
      boundary = 'bottom right top'
    []

    [inject_heat]
        type = DirichletBC
        variable = temperature
        # function = ${inject_temp}
        boundary = 'injection_area'
        value = 45 #318.15
    []

    [produce_heat]
        type = PorousFlowSink
        variable = temperature
        boundary = injection_area
        flux_function = production_rate_value
        fluid_phase = 0
        use_enthalpy = true
        #save_in = heat_flux_out
    []

    [inject_fluid]
      type = PorousFlowSink
      variable = porepressure
      boundary = injection_area
      flux_function = injection_rate_value
    []

    [produce_fluid]
      type = PorousFlowSink
      variable = porepressure
      boundary = injection_area
      flux_function = production_rate_value
    []
[]

[Controls]
    [inject_on]
      type = ConditionalFunctionEnableControl
      enable_objects = 'BCs::inject_fluid'
      conditional_function = inject
      implicit = false
      execute_on = 'initial timestep_begin'
    []
    [produce_on]
      type = ConditionalFunctionEnableControl
      enable_objects = 'BCs::produce_fluid'
      conditional_function = produce
      implicit = false
      execute_on = 'initial timestep_begin'
    []

    [temp_on]
      type = ConditionalFunctionEnableControl
      enable_objects = 'BCs::inject_heat'
      conditional_function = inject
      execute_on = 'initial timestep_begin'
      implicit = false
    []
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
      symbol_values = true_screen_area
      expression = '-${inject_fluid_mass}/(true_screen_area * ${inject_time})'
    []
    [production_rate_value]
      type = ParsedFunction
      symbol_names = true_screen_area
      symbol_values = true_screen_area
      expression = '${produce_fluid_mass}/(true_screen_area * ${produce_time})'
    []
  
[]


[FluidProperties]
    [true_water]
      type = HydrogenFluidProperties # SimpleFluidProperties  # #
    #   bulk_modulus = 14.53e6
    #   density0 = 9.55 # 905.24
    #   viscosity = 9.93e-6
    #   cv = 14657.2
    #   cp = 14657.2
    #   porepressure_coefficient = 0.0
    #   thermal_expansion = 0
    #   thermal_conductivity = 0.22
    []
    #[hydrogen]
    #  type = HydrogenFluidProperties # BrineFluidProperties #   HydrogenFluidProperties # BrineFluidProperties # #HydrogenFluidProperties #BrineFluidProperties
    #[]

[]

[Materials]

    [temperature]
      type = PorousFlowTemperature
      temperature = temperature
      #at_nodes = true
    []
    [ppss]
        type = PorousFlow1PhaseFullySaturated
        porepressure = porepressure
    []

    [massfrac]
        type = PorousFlowMassFraction
        mass_fraction_vars = ''
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
      block = aquifer
    []

    [porosity_caps]
        type = PorousFlowPorosityConst
        block = caps
        porosity = 0.1
    []

    [permeability_reservoir]
      type = PorousFlowPermeabilityConst
      #permeability = '1.8e-11 0 0  0 1.8e-11 0  0 0 1.8e-11' #3e-13 = 300 mD
      permeability = '6.6646e-13 0 0  0 6.6646e-13 0  0 0 6.6646e-13'
      block = aquifer
    []

    [permeability_caps]
        type = PorousFlowPermeabilityConst
        block = caps
        permeability = '5E-18 0 0   0 5E-18 0   0 0 0'
    [] 

    [relperm]
     type = PorousFlowRelativePermeabilityConst
     phase = 0
     #at_nodes = true
    []

    [thermal_conductivity_reservoir]
      type = PorousFlowThermalConductivityIdeal
     
      dry_thermal_conductivity = '3.083 0 0  0 3.083 0  0 0 0'
      #wet_thermal_conductivity = '3.083 0 0  0 0 0  0 0 3.083'
      #wet_thermal_conductivity = '0 0 0  0 0 0  0 0 3.083'
      block = aquifer
    []

    [thermal_conductivity_caps]
        type = PorousFlowThermalConductivityIdeal
        block = caps
        dry_thermal_conductivity = '4.310 0 0  0 4.310 0  0 0 0'
        #wet_thermal_conductivity = '0 0 0  0 4.310 0  0 0 0'
    []

    [internal_energy_reservoir]
      type = PorousFlowMatrixInternalEnergy
      specific_heat_capacity = 951.6
      density = 2601.2
      block = aquifer
    []

    [internal_energy_caps]
        type = PorousFlowMatrixInternalEnergy
        block = caps
        specific_heat_capacity = 829
        density = 2773.4
    []

[]

[UserObjects]
    [dictator]
      type = PorousFlowDictator
      porous_flow_vars = 'porepressure temperature'
      number_fluid_phases = 1
      number_fluid_components = 1
    []
    [pc]
      type = PorousFlowCapillaryPressureConst
      pc = 0
    []

    [fluid_advective_flux]
      type = PorousFlowAdvectiveFluxCalculatorSaturated
      flux_limiter_type = superbee
    []
    [heat_advective_flux]
      type = PorousFlowAdvectiveFluxCalculatorSaturatedHeat
      flux_limiter_type = superbee
    []


[]
  
[Kernels]
    [mass_gas_dot]
      type = PorousFlowMassTimeDerivative
      fluid_component = 0
      variable = porepressure
      multiply_by_density = true
      
    []
    [flux_gas]
      type =  PorousFlowFluxLimitedTVDAdvection #PorousFlowAdvectiveFlux # PorousFlowFullySaturatedDarcyFlow
      #fluid_component = 0
      use_displaced_mesh = false
      variable = porepressure
      #gravity = '0 0 0'
      advective_flux_calculator = fluid_advective_flux
    []
    [energy_dot]
      type = PorousFlowEnergyTimeDerivative
      variable = temperature
    []
    [advection]
      type = PorousFlowFluxLimitedTVDAdvection #PorousFlowHeatAdvection
      use_displaced_mesh = false
      variable = temperature
      advective_flux_calculator = heat_advective_flux
    []
    [conduction]
      type = PorousFlowHeatConduction #PorousFlowExponentialDecay
      use_displaced_mesh = false
      variable = temperature
      #reference = 331
      #rate = rate
    []
[]
  
  
[Postprocessors]

    [true_screen_area] # this accounts for meshes that do not match screen_top and screen_bottom exactly
    type = AreaPostprocessor
    boundary = injection_area
    execute_on = 'initial'
    outputs = 'none'
    []

  
[]
  
[VectorPostprocessors]
    [ptsussfull_res1]
      type = LineValueSampler
      use_displaced_mesh = false
      start_point = '0.1 5 0'
      end_point = '100 5 0'
      sort_by = x
      num_points = 1000
      outputs = csv
      variable = 'porepressure temperature'
    []

    [ptsussfull_res2]
        type = LineValueSampler
        use_displaced_mesh = false
        start_point = '0.1 0 0'
        end_point = '100 0 0'
        sort_by = x
        num_points = 1000
        outputs = csv
        variable = 'porepressure temperature'
    []


    [ptsussfull_cap1]
            type = LineValueSampler
            use_displaced_mesh = false
            start_point = '0.1 10 0' #5m above res
            end_point = '100 10 0'
            sort_by = x
            num_points = 1000
            outputs = csv
            variable = 'porepressure temperature'
    []

    [ptsussfull_cap2]
        type = LineValueSampler
        use_displaced_mesh = false
        start_point = '0.1 15 0' #5m above res
        end_point = '100 15 0'
        sort_by = x
        num_points = 1000
        outputs = csv
        variable = 'porepressure temperature'
    []

    [ptsussfull_cap3]
        type = LineValueSampler
        use_displaced_mesh = false
        start_point = '0.1 20 0' #5m above res
        end_point = '100 20 0'
        sort_by = x
        num_points = 1000
        outputs = csv
        variable = 'porepressure temperature'
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
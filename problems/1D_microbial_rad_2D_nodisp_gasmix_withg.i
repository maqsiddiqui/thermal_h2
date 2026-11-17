inject_fluid_mass = 3614484 # kg
inject_fluid_mass_h = '${fparse inject_fluid_mass*0.13/1.02}' # kg
inject_fluid_mass_c = '${fparse inject_fluid_mass*0.87/1.02}' # kg
inject_time = 1.577e7 # # 6 months in seconds
store_time =  0 # seconds
produce_time = 0 #1.577e7 # # 6 months in seconds
rest_time = 0 # days
num_cycles = 1 # Currently needs to be <= 10

cycle_length = '${fparse inject_time + store_time + produce_time + rest_time}'

end_simulation = '${fparse cycle_length * num_cycles}'

start_injection1 = 0
three_months = 7.884e6 #secs
one_month = 2.628e6 #secs
end_injection1 = '${fparse start_injection1 + inject_time}'

synctimes = '${start_injection1} ${one_month} ${three_months} ${end_injection1}'


[Mesh]
    construct_side_list_from_node_list=true
    [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1000
    # bias_x = 1.05
    xmin = 0.1
    xmax = 1000
    ny = 5
    ymin = 0
    ymax = 5
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
    gravity = '0 -9.81 0'
[]

[UserObjects]
    [dictator]
      type = PorousFlowDictator
      porous_flow_vars = 'porepressure cgH2 cgCH4 cgCO2'
      number_fluid_phases = 1
      number_fluid_components = 4
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
      initial_condition = 0#2170 psi
      scaling = 1e4
    []
    [cgCH4]
      initial_condition = 0 #
      scaling = 1e4
    []
    [cgCO2]
      initial_condition = 0 #
      scaling = 1e4
    []
    [npop_d]
      initial_condition = 1 #
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

    [source_h2]
        order = FIRST
        family = LAGRANGE
    []

    [source_ch4]
        order = FIRST
        family = LAGRANGE
    []

    [source_co2]
        order = FIRST
        family = LAGRANGE
    []

  [density]
    family = MONOMIAL
    order = FIRST
  []
[]

[Kernels]

    [npop_dt]
        type = TimeDerivative
        variable = npop_d
    []

    [npop_balance]
        type = MatReaction
        variable = npop_d
        reaction_rate = 'npop_source'
        
    []

    [mass_gas_dot]
        type = PorousFlowMassTimeDerivative
        fluid_component = 3
        variable = porepressure
        multiply_by_density = true
    []

    [flux_gas]
        type = PorousFlowFullySaturatedDarcyFlow
        fluid_component = 3
        variable = porepressure
        multiply_by_density = true
    []

    # [diff_gas]
    #     type = PorousFlowDispersiveFlux
    #     fluid_component = 0
    #     variable = porepressure
    #     disp_trans = 0
    #     disp_long = 0
    # []

    [cgH2_dt]
        type = PorousFlowMassTimeDerivative
        fluid_component = 0
        variable = cgH2
        multiply_by_density = true
    []
    # [cgH2_diff]
    #     type = PorousFlowDispersiveFlux
    #     fluid_component = 1
    #     variable = cgH2
    #     disp_long = 0
    #     disp_trans = 0
    # []
    [cgH2_conv]
        type =  PorousFlowFullySaturatedDarcyFlow
        fluid_component = 0
        variable = cgH2
        multiply_by_density = true
    []

    [cgH2_reaction]
        type = MatReaction
        variable = cgH2
        reaction_rate = 'doublemonodterm_h2'
        
    []

    [cgCH4_dt]
        type = PorousFlowMassTimeDerivative
        fluid_component = 1
        variable = cgCH4
        multiply_by_density = true
    []
    # [cgCH4_diff]
    #     type = PorousFlowDispersiveFlux
    #     fluid_component = 2
    #     variable = cgCH4
    #     disp_long = 0
    #     disp_trans = 0

    # []
    [cgCH4_conv]
        type = PorousFlowFullySaturatedDarcyFlow
        fluid_component = 1
        variable = cgCH4
        multiply_by_density = true
    []

    [cgCH4_reaction]
        type = MatCoupledForce
        variable = cgCH4
        coef = 1 #
        v = cgH2
        material_properties = 'doublemonodterm_ch4'
        
    []

    [cgCO2_dt]
        type = PorousFlowMassTimeDerivative
        fluid_component = 2
        variable = cgCO2
        multiply_by_density = true
    []
    # [cgCO2_diff]
    #     type = PorousFlowDispersiveFlux
    #     fluid_component = 3
    #     variable = cgCO2
    #     disp_long = 0
    #     disp_trans = 0
    # []
    [cgCO2_conv]
        type = PorousFlowFullySaturatedDarcyFlow
        fluid_component = 2
        variable = cgCO2
        multiply_by_density = true
    []

    [cgCO2_reaction]
        type = MatCoupledForce
        variable = cgCO2
        coef = -1 #
        v = cgH2
        material_properties = 'doublemonodterm_co2'
        
    []

[]
  
[AuxKernels]
    # [source_h2]
    #     type = ParsedAux
    #     variable = source_h2
    #     coupled_variables = 'cgH2'
    #     expression = '1.59e-8*cgH2' # 
    # []

    # [source_ch4]
    #     type = ParsedAux
    #     variable = source_ch4
    #     coupled_variables = 'cgH2'
    #     expression = '3.18e-8*cgH2'
    # []

    [density]
      type = PorousFlowPropertyAux
      variable = density
      property = density
      phase = 0
      execute_on = 'initial timestep_end'
    []

[]
  

[Functions]
    [injection_rate]
      type = ParsedFunction
      symbol_values = injection_area
      symbol_names = area
      expression = '-${inject_fluid_mass}/(area *0.162 * ${inject_time})'
    []

    [H2_injection_rate]
      type = ParsedFunction
      symbol_values = injection_area
      symbol_names = area
    #   expression = '-${inject_fluid_mass}*0.13/(area *0.162*1.02 * ${inject_time})' #Retard_h = 1.02
      expression = '-${inject_fluid_mass_h}/(area *0.162 * ${inject_time})' #Retard_h = 1.02
    []
    [CO2_injection_rate]
      type = ParsedFunction
      symbol_values = injection_area
      symbol_names = area
    #   expression = '-${inject_fluid_mass}*0.87/(area *0.162 *3.05 * ${inject_time})'
      expression = '-${inject_fluid_mass_c}/(area *0.162 * ${inject_time})'
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
    # [water]
    #  type = HydrogenFluidProperties # BrineFluidProperties #   HydrogenFluidProperties # BrineFluidProperties # #HydrogenFluidProperties #BrineFluidProperties
    # []

    [CH4]
      type = MethaneFluidProperties # SimpleFluidProperties  # #
    []
    [H2]
      type = HydrogenFluidProperties # BrineFluidProperties #   HydrogenFluidProperties # BrineFluidProperties # #HydrogenFluidProperties #BrineFluidProperties
    []
    [N2]
    type = NitrogenFluidProperties
    []
    [CO2]
    type = CO2FluidProperties
    []

[]
  
[Materials]

    [npop_source]
        type = ParsedMaterial
        property_name = npop_source
        constant_names = 'n0 Pgrowth_max Y_h2 Sw Mg Mratio poros gamma halfsat_h2 halfsat_co2 Hc_h Hc_c Ret decay_rate'
        constant_expressions = '1.8e14 1.95e-4 4.2e11 0.1 0.002 0.022 0.18 1 2e-8 4.4e-5 0.02 1.59 1.02 9e-8'
        coupled_variables = 'cgH2 cgCO2 npop_d'
        expression = '(npop_d*Pgrowth_max*Hc_c*Hc_h*cgCO2*cgH2/Y_h2/(halfsat_h2+(cgH2*Hc_h))/(halfsat_co2+(cgCO2*Hc_c))) - (decay_rate*npop_d)'
        #outputs = exodus
    []

    [doublemonodsource_h2]
        type = ParsedMaterial
        property_name = doublemonodterm_h2
        constant_names = 'n0 Pgrowth_max Y_h2 Sw Mg Mratio poros gamma halfsat_h2 halfsat_co2 Hc_h Hc_c Ret'
        constant_expressions = '1.8e14 1.95e-4 4.2e11 0.1 0.002 0.022 0.18 1 2e-8 4.4e-5 0.02 1.59 1.02'
        coupled_variables = 'cgH2 cgCO2 npop_d'
        expression = '-gamma*Mratio*poros*Sw*n0*npop_d*Pgrowth_max*Mg*Hc_c*Hc_h*cgCO2/Y_h2/Ret/(halfsat_h2+(cgH2*Hc_h))/(halfsat_co2+(cgCO2*Hc_c))'
        #outputs = exodus
    []

    [doublemonodsource_co2]
        type = ParsedMaterial
        property_name = doublemonodterm_co2
        constant_names = 'n0 Pgrowth_max Y_h2 Sw Mg Mratio poros gamma halfsat_h2 halfsat_co2 Hc_h Hc_c Ret'
        constant_expressions = '1.8e14 1.95e-4 4.2e11 0.1 0.002 0.022 0.18 5.5 2e-8 4.4e-5 0.02 1.59 3.05'
        coupled_variables = 'cgH2 cgCO2 npop_d'
        expression = 'gamma*Mratio*poros*Sw*n0*npop_d*Pgrowth_max*Mg*Hc_c*Hc_h*cgCO2/Y_h2/Ret/(halfsat_h2+(cgH2*Hc_h))/(halfsat_co2+(cgCO2*Hc_c))'
        #outputs = exodus
    []

    [doublemonodsource_ch4]
        type = ParsedMaterial
        property_name = doublemonodterm_ch4
        constant_names = 'n0 Pgrowth_max Y_h2 Sw Mg Mratio poros gamma halfsat_h2 halfsat_co2 Hc_h Hc_c Ret'
        constant_expressions = '1.8e14 1.95e-4 4.2e11 0.1 0.002 0.022 0.18 2 2e-8 4.4e-5 0.02 1.59 1.06'
        coupled_variables = 'cgH2 cgCO2 npop_d'
        expression = 'gamma*Mratio*poros*Sw*n0*npop_d*Pgrowth_max*Mg*Hc_c*Hc_h*cgCO2/Y_h2/Ret/(halfsat_h2+(cgH2*Hc_h))/(halfsat_co2+(cgCO2*Hc_c))'
        #outputs = exodus
    []

    [temperature]
      type = PorousFlowTemperature
      temperature = 341.15
      #at_nodes = true
    []
    [ppss]
        type = PorousFlow1PhaseFullySaturated
        porepressure = porepressure
    []

    [H2_CH4mix]
        type = PorousFlowMultiComponentGasMixture
        fp = 'H2 CH4 CO2 N2'
        mass_fraction = 'cgH2 cgCH4 cgCO2'
        phase = 0
        #compute_enthalpy = true
        #compute_internal_energy = true
        #temperature_unit = Celsius
    []

    [massfrac]
        type = PorousFlowMassFraction
        mass_fraction_vars = 'cgH2 cgCH4 cgCO2'
    []

    # [water]
    #   type =  PorousFlowSingleComponentFluid # PorousFlowMultiComponentFluid # PorousFlowSingleComponentFluid #
    #   fp = water
    #   #fp = H2
    #   #x = xH2
    #   phase = 0
    #   #at_nodes = true
    # []

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
     s_res = 0.1
     sum_s_res = 0.1
     #at_nodes = true
    []

    # [diff]
    #     type = PorousFlowDiffusivityConst
    #     diffusion_coeff = '4.4e-5 4.4e-5 4.4e-5'
    #     tortuosity = 1
    # []

  []
  
[BCs]

    # [cgH2_in]
    #     type = DirichletBC
    #     variable = cgH2
    #     boundary = injection_area
    #     value = 0.13
    # []
    # [cgCO2_in]
    #       type = DirichletBC
    #       variable = cgCO2
    #       boundary = injection_area
    #       value = 0.87
    # []
    # [cgCH4_in]
    #       type = DirichletBC
    #       variable = cgCH4
    #       boundary = injection_area
    #       value = 0
    # []


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

    [total_injection]
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

    [H2_injection]
      type = PorousFlowSink
      boundary = injection_area
      variable = cgH2
      use_mobility = false
      use_relperm = false
      flux_function =  H2_injection_rate
      #value = -0.000125
    #   mass_fraction_component = 0
    #   fluid_phase=0
    []

    [CO2_injection]
      type = PorousFlowSink
      boundary = injection_area
      variable = cgCO2
      use_mobility = false
      use_relperm = false
      flux_function =  CO2_injection_rate
      #value = -0.000125
    #   mass_fraction_component = 2
    #   fluid_phase=0
    []

    # [inject_right]
    #     type = PorousFlowOutflowBC
    #     variable = porepressure
    #     boundary = right
    #     mass_fraction_component = 0
    #     # flux_function = injection_rate_value
    # []

    #[outflow_a]
    #    type = PorousFlowOutflowBC
    #    boundary = right
    #    include_relperm = false # no need for relperm in this fully-saturated simulation
    #    #mass_fraction_component = 0
    #   variable = pgas
    #  []

[]
  
[Postprocessors]

    [injection_area]
        type = AreaPostprocessor
        boundary = injection_area
        execute_on = initial
    []

    [injection_flux]
        type = FunctionValuePostprocessor
        function = injection_rate
        execute_on = initial
    []
    

[]
  
[VectorPostprocessors]
    [ptsuss]
      type = LineValueSampler
      use_displaced_mesh = false
      start_point = '0.1 2.5 0'
      end_point = '1000 2.5 0'
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


# [Preconditioning]
#   [basic]
#     type = SMP
#     full = true
#     petsc_options_iname = '-pc_type -pc_factor_shift_type '
#     petsc_options_value = ' lu       NONZERO'
#   []
# []
  
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
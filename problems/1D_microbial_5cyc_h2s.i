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

[Variables]
    [cgH2]
      order = FIRST
      family = LAGRANGE
    []
    [cgCH4]
        order = FIRST
        family = LAGRANGE
    []
    [cgCO2]
        order = FIRST
        family = LAGRANGE
    []
    [cgH2S]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxVariables]
    [pressure]
      order = FIRST
      family = LAGRANGE
    []
    [source_h2]
        order = FIRST
        family = LAGRANGE
    []    
    [source_ch4]
        order = FIRST
        family = LAGRANGE
    []
    [source_CO2]
        order = FIRST
        family = LAGRANGE
    []
    [source_H2S]
        order = FIRST
        family = LAGRANGE
    []
[]

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
    [cgCO2]
        type = ConstantIC
        variable = cgCO2
        value = 0
    []
    [cgH2S]
        type = ConstantIC
        variable = cgH2S
        value = 0
    []
    [pressure]
        type = FunctionIC
        variable = pressure
        function = 14.55e6-x
    []
[]

[Kernels]
    [cgH2_dt]
      type = TimeDerivative
      variable = cgH2
    []
    [cgH2_diff]
        type = MatDiffusion
        variable = cgH2
        diffusivity = 4.4e-5
    []
    [cgH2_conv]
        type = ConservativeAdvection
        variable = cgH2
        upwinding_type = full
        velocity = '6.2e-6 0 0'
    []

    [cgH2_ex]
      type = ConservativeAdvection #PrimaryConvection
      variable = cgH2
      upwinding_type = full
      velocity = '-6.2e-6 0 0'
      #p = pressure
    []

    [cgH2_reaction]
        type = CoupledForce
        coef = -1 #
        v = source_h2
        variable = cgH2
        
    []

    [cgCO2_dt]
        type = TimeDerivative
        variable = cgCO2
    []
    [cgCO2_diff]
          type = MatDiffusion
          variable = cgCO2
          diffusivity = 4.4e-5
    []
    [cgCO2_conv]
          type = ConservativeAdvection
          variable = cgCO2
          upwinding_type = full
          velocity = '6.2e-6 0 0'
    []

    [cgCO2_reaction]
        type = CoupledForce
        coef = -1 #
        v = source_CO2
        variable = cgCO2
        
    []

    [cgCO2_ex]
        type = ConservativeAdvection #PrimaryConvection
        variable = cgCO2
        upwinding_type = full
        velocity = '-6.2e-6 0 0'
        #p = pressure
    []

    [cgCH4_dt]
        type = TimeDerivative
        variable = cgCH4
    []
    [cgCH4_diff]
          type = MatDiffusion
          variable = cgCH4
          diffusivity = 4.4e-5
    []
    [cgCH4_conv]
      type =  ConservativeAdvection #PrimaryConvection
      variable = cgCH4
      velocity = '6.2e-6 0 0'
      upwinding_type = full
      #p = pressure
    []
    [cgCH4_ex]
      type =  ConservativeAdvection #PrimaryConvection
      variable = cgCH4
      velocity = '-6.2e-6 0 0'
      upwinding_type = full
      #p = pressure
    []

    [cgCH4_reaction]
        type = CoupledForce
        coef = 1 #
        v = source_ch4
        variable = cgCH4
    []

    [cgH2S_dt]
        type = TimeDerivative
        variable = cgH2S
    []

    [cgH2S_diff]
        type = MatDiffusion
        variable = cgH2S
        diffusivity = 4.4e-5
    []
    [cgH2S_conv]
        type =  ConservativeAdvection #PrimaryConvection
        variable = cgH2S
        velocity = '6.2e-6 0 0'
        upwinding_type = full
        #p = pressure
    []
    [cgH2S_ex]
        type =  ConservativeAdvection #PrimaryConvection
        variable = cgH2S
        velocity = '-6.2e-6 0 0'
        upwinding_type = full
        #p = pressure
    []

    [cgH2S_reaction]
        type = CoupledForce
        coef = 1 #
        v = source_H2S
        variable = cgH2S
    []

[]

[AuxKernels]

    [source_h2]
        type = ParsedAux
        variable = source_h2
        coupled_variables = 'cgH2 cgCO2'
        expression = '9.34e-8*cgCO2*cgH2 + 6.16e-9*cgH2'
    []
    [source_CO2]
        type = ParsedAux
        variable = source_CO2
        coupled_variables = 'cgH2 cgCO2'
        expression = '0.25*9.34e-8*cgCO2*cgH2'
    []
    [source_ch4]
        type = ParsedAux
        variable = source_ch4
        coupled_variables = 'cgH2 cgCO2'
        expression = '0.25*9.34e-8*cgCO2*cgH2'
    []
    [source_H2S]
        type = ParsedAux
        variable = source_H2S
        coupled_variables = 'cgH2'
        expression = '0.2*6.16e-9*cgH2'
    []
[]

[BCs]
    [cgH2_left]
      type = DirichletBC
      variable = cgH2
      boundary = left
      value = 0.8
    []

    [cgCO2_left]
        type = DirichletBC
        variable = cgCO2
        boundary = left
        value = 0.2
    []

    [cgCH4_left]
        type = DirichletBC
        variable = cgCH4
        boundary = left
        value = 0
    []

    [cgH2S_left]
        type = DirichletBC
        variable = cgH2S
        boundary = left
        value = 0
    []
    

[]

[Controls]
    [inject_h2_bc]
      type = ConditionalFunctionEnableControl
      enable_objects = 'BCs::cgH2_left'
      conditional_function = inject
      implicit = false
      execute_on = 'initial timestep_begin'
    []

    [prod_h2_kern]
      type = ConditionalFunctionEnableControl
      enable_objects = 'Kernels::cgH2_ex'
      conditional_function = produce
      implicit = false
      execute_on = 'initial timestep_begin'
    []

    [inject_ch4_bc]
      type = ConditionalFunctionEnableControl
      enable_objects = 'BCs::cgCH4_left'
      conditional_function = inject
      implicit = false
      execute_on = 'initial timestep_begin'
    []

    [prod_ch4_kern]
      type = ConditionalFunctionEnableControl
      enable_objects = 'Kernels::cgCH4_ex'
      conditional_function = produce
      implicit = false
      execute_on = 'initial timestep_begin'
    []

    [inject_co2_bc]
        type = ConditionalFunctionEnableControl
        enable_objects = 'BCs::cgCO2_left'
        conditional_function = inject
        implicit = false
        execute_on = 'initial timestep_begin'
    []
  
    [prod_co2_kern]
        type = ConditionalFunctionEnableControl
        enable_objects = 'Kernels::cgCO2_ex'
        conditional_function = produce
        implicit = false
        execute_on = 'initial timestep_begin'
    []

    [inject_h2s_bc]
        type = ConditionalFunctionEnableControl
        enable_objects = 'BCs::cgH2S_left'
        conditional_function = inject
        implicit = false
        execute_on = 'initial timestep_begin'
    []
  
    [prod_h2s_kern]
        type = ConditionalFunctionEnableControl
        enable_objects = 'Kernels::cgH2S_ex'
        conditional_function = produce
        implicit = false
        execute_on = 'initial timestep_begin'
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

[]

# [Materials]
#     [porous_inj]
#       type = GenericConstantMaterial
#       prop_names = 'diffusivity conductivity porosity'
#       prop_values = '4.4e-5 3.27e-5 1'
#     []
#     [porous_prod]
#         type = GenericConstantMaterial
#         prop_names = 'diffusivityp conductivityp porosityp'
#         prop_values = '4.4e-5 -3.27e-5 1'
#     []
# []

[VectorPostprocessors]
    [lineplot]
      type = LineValueSampler
      use_displaced_mesh = false
      start_point = '0 0 0'
      end_point = '1000 0 0'
      sort_by = x
      num_points = 200
      outputs = csv
      variable = 'cgH2 cgCH4 cgCO2 cgH2S pressure'
    []

    [wellpoint]
      type = PointValueSampler
      use_displaced_mesh = false
      points = '0 0 0'
      sort_by = x
      outputs = csv
      variable = 'cgH2 cgCH4 cgCO2 cgH2S pressure'
      execute_on = 'timestep_begin'
      # contains_complete_history = true
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
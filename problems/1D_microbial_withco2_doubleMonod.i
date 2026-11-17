inject_time = 4.32e7 # # 500 days in seconds
store_time =  0 # seconds
produce_time = 0 # # 6 months in seconds
rest_time = 0 # days
num_cycles = 1 # Currently needs to be <= 10

cycle_length = '${fparse inject_time + store_time + produce_time + rest_time}'

end_simulation = '${fparse cycle_length * num_cycles}'

start_injection1 = 0
three_months = 7.884e6 #secs
one_month = 2.628e6 #secs
end_injection1 = '${fparse start_injection1 + inject_time}'

synctimes = '${start_injection1} ${one_month} ${three_months} ${end_injection1}'


#end_simulation = 1.577e7
#synctimes = '2.628e6 7.884e6 1.577e7' # 1 mo, 3 mo, 6 mo

[Mesh]
    type = GeneratedMesh
    dim = 1
    xmin = 0
    xmax = 1000
    nx = 200
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
    # [cgH2_diff]
    #     type = MatDiffusion
    #     variable = cgH2
    #     diffusivity = 4.4e-5
    # []
    [cgH2_conv]
        type = ConservativeAdvection
        variable = cgH2
        upwinding_type = full
        velocity = '3.27e-5 0 0'
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
    # [cgCO2_diff]
    #       type = MatDiffusion
    #       variable = cgCO2
    #       diffusivity = 4.4e-5
    # []
    [cgCO2_conv]
          type = ConservativeAdvection
          variable = cgCO2
          upwinding_type = full
          velocity = '3.27e-5 0 0'
    []

    [cgCO2_reaction]
        type = CoupledForce
        coef = -1 #
        v = source_CO2
        variable = cgCO2
        
    []

    [cgCH4_dt]
        type = TimeDerivative
        variable = cgCH4
    []
    # [cgCH4_diff]
    #       type = MatDiffusion
    #       variable = cgCH4
    #       diffusivity = 4.4e-5
    # []
    [cgCH4_conv]
      type =  ConservativeAdvection #PrimaryConvection
      variable = cgCH4
      velocity = '3.27e-5 0 0'
      upwinding_type = full
      #p = pressure
    []


    [cgCH4_reaction]
        type = CoupledForce
        coef = 1 #
        v = source_ch4
        variable = cgCH4
    []

[]

[AuxKernels]

    [source_h2]
        type = ParsedAux
        variable = source_h2
        coupled_variables = 'cgH2 cgCO2'
        expression = '9.34e-8*(cgCO2/(0.1+cgCO2))*(cgH2/(0.1+cgH2))'
    []
    [source_CO2]
        type = ParsedAux
        variable = source_CO2
        coupled_variables = 'cgH2 cgCO2'
        expression = '0.25*9.34e-8*(cgCO2/(0.1+cgCO2))*(cgH2/(0.1+cgH2))'
    []
    [source_ch4]
        type = ParsedAux
        variable = source_ch4
        coupled_variables = 'cgH2 cgCO2'
        expression = '0.25*9.34e-8*(cgCO2/(0.1+cgCO2))*(cgH2/(0.1+cgH2))'
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

[]

# [Controls]

# []

[Functions]



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
      end_point = '500 0 0'
      sort_by = x
      num_points = 100
      outputs = csv
      variable = 'cgH2 cgCH4 cgCO2 pressure'
    []

    # [wellpoint]
    #   type = PointValueSampler
    #   use_displaced_mesh = false
    #   points = '0 0 0'
    #   sort_by = x
    #   outputs = csv
    #   variable = 'cgH2 cgCH4 cgCO2 pressure'
    #   execute_on = 'timestep_begin'
    #   # contains_complete_history = true
    # []


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
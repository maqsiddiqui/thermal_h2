end_simulation = 1.577e7
synctimes = '2.628e6 7.884e6 1.577e7' # 1 mo, 3 mo, 6 mo

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
[]

[AuxVariables]
    [pressure]
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
    [pressure]
        type = FunctionIC
        variable = pressure
        function = 5e6-x
    []
[]

[Kernels]
    [cgH2_dt]
      type = PrimaryTimeDerivative
      variable = cgH2
    []
    [cgH2_diff]
        type = PrimaryDiffusion
        variable = cgH2
    []
    [cgH2_conv]
        type = PrimaryConvection
        variable = cgH2
        p = pressure
    []

    [cgH2_reaction]
        type = CoefReaction
        coefficient = 9.34e-8
        variable = cgH2
    []

    [cgCH4_dt]
        type = PrimaryTimeDerivative
        variable = cgCH4
    []
    [cgCH4_diff]
          type = PrimaryDiffusion
          variable = cgCH4
    []
    [cgCH4_conv]
          type = PrimaryConvection
          variable = cgCH4
          p = pressure
    []

    [cgCH4_reaction]
        type = CoupledForce
        coef = 2.34e-8 #
        v = cgH2
        variable = cgCH4

    []

[]

[BCs]
    [cgH2_left]
      type = DirichletBC
      variable = cgH2
      boundary = left
      value = 0.8
    []

    [cgCH4_left]
        type = DirichletBC
        variable = cgCH4
        boundary = left
        value = 0
    []

[]

[Materials]
    [porous]
      type = GenericConstantMaterial
      prop_names = 'diffusivity conductivity porosity'
      prop_values = '4.19e-5 3.1e-5 1'
    []
[]

[VectorPostprocessors]
    [1D_microbial]
      type = LineValueSampler
      use_displaced_mesh = false
      start_point = '0 0 0'
      end_point = '1000 0 0'
      sort_by = x
      num_points = 50
      outputs = csv
      variable = 'cgH2 cgCH4 pressure'
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
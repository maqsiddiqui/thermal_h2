time = 1.577e7

synctimes = '2.628e6 7.884e6 1.577e7'


aq_thickness = 50
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
        dx = '1 2 12 7 8'
        ix = '10 10 50 10 10'
        dy = '25 70 25'
        iy = '10 70 10'
    []

    [shift]
        type = TransformGenerator
        input = mesh
        transform = translate
        vector_value = '0.1 -60 0'
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
    gravity = '0 0 0'
[]

[Variables]
    [porepressure]
    []
    [mass_frac_CH4]
        scaling = 1e4
    []
[]

[AuxVariables]
    [density]
      family = MONOMIAL
      order = FIRST
    []
[]


[ICs]

    [mass_frac_CH4]
        type = FunctionIC
        variable = mass_frac_CH4
        function = s_curve
        block = aquifer
    []
    [pressure]
      type = ConstantIC
      variable = porepressure
      value = 5e6
    []

[]

[Functions]
  [s_curve]
      type = ParsedFunction
      # expression = 'if(y > (x <= 30 ? 12.5 + 12.5*cos(pi*(1 - x/30)):25), 0, 1)'
      expression = 'if(y > 12.5 + 12.5*cos(pi*(1 - x/30)), 0, 1)'
      # vars = 'k m'
      # vals = '5 5'
  []

[]


# [BCs]

#     [outer_boundary_porepressure]
#       type = DirichletBC
#       preset = true
#       variable = porepressure
#     #   function = insitu_pressure
#       value = 5e6
#       boundary = 'left bottom right top'
#     []

# []

# [Controls]
#     [inject_on]
#       type = ConditionalFunctionEnableControl
#       enable_objects = 'BCs::inject_fluid'
#       conditional_function = inject
#       implicit = false
#       execute_on = 'initial timestep_begin'
#     []
#     [produce_on]
#       type = ConditionalFunctionEnableControl
#       enable_objects = 'BCs::produce_fluid'
#       conditional_function = produce
#       implicit = false
#       execute_on = 'initial timestep_begin'
#     []

# []


[FluidProperties]
    [CH4]
      type = MethaneFluidProperties # SimpleFluidProperties  # #
    []
    [H2]
      type = HydrogenFluidProperties # BrineFluidProperties #   HydrogenFluidProperties # BrineFluidProperties # #HydrogenFluidProperties #BrineFluidProperties
    []

[]

[Materials]

    [temperature]
      type = PorousFlowTemperature
      temperature = 313
      #at_nodes = true
    []
    [ppss]
        type = PorousFlow1PhaseFullySaturated
        porepressure = porepressure
    []

    [H2]
        type = PorousFlowMultiComponentGasMixture
        fp = 'CH4 H2'
        mass_fraction = mass_frac_CH4
        phase = 0
        compute_enthalpy = false
        compute_internal_energy = false
        # temperature_unit = Celsius
    []

    [massfrac]
        type = PorousFlowMassFraction
        mass_fraction_vars = mass_frac_CH4
    []

    [diff]
        type = PorousFlowDiffusivityConst
        diffusion_coeff = '1e-7 1e-7'
        tortuosity = 1
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

[]

[UserObjects]
    [dictator]
      type = PorousFlowDictator
      porous_flow_vars = 'porepressure mass_frac_CH4'
      number_fluid_phases = 1
      number_fluid_components = 2
    []
    [pc]
      type = PorousFlowCapillaryPressureConst
      pc = 0
    []

[]
  
[Kernels]
    [mass_gas_dot]
      type = PorousFlowMassTimeDerivative
      fluid_component = 1
      variable = porepressure
      multiply_by_density = true
      
    []
    [flux_gas]
      type =  PorousFlowFullySaturatedDarcyFlow #PorousFlowFluxLimitedTVDAdvection #PorousFlowAdvectiveFlux # PorousFlowFullySaturatedDarcyFlow
      fluid_component = 1
      use_displaced_mesh = false
      variable = porepressure
      #gravity = '0 0 0'
      # advective_flux_calculator = fluid_advective_flux

    []
    [diff_CO2]
        type = PorousFlowDispersiveFlux
        fluid_component = 1
        variable = porepressure
        disp_trans = 0
        disp_long = 0
    []
    
    [mass_CH4]
        type = PorousFlowMassTimeDerivative
        fluid_component = 0
        variable = mass_frac_CH4
    []
    [adv_CH4]
        type = PorousFlowFullySaturatedDarcyFlow
        fluid_component = 0
        variable = mass_frac_CH4
    []
    [diff_CH4]
        type = PorousFlowDispersiveFlux
        fluid_component = 0
        variable = mass_frac_CH4
        disp_trans = 0
        disp_long = 0
    []

[]
  
[AuxKernels]
    [density]
      type = PorousFlowPropertyAux
      variable = density
      property = density
      phase = 0
      execute_on = 'initial timestep_end'
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
    [mass_frac_CH4]
        type = LineValueSampler
        start_point = '0.1 25 0'
        end_point = '30 25 0'
        variable = mass_frac_CH4
        sort_by = x
        num_points = 300
    []

    # [ptsussfull_res1]
    #   type = LineValueSampler
    #   use_displaced_mesh = false
    #   start_point = '0.1 25 0'
    #   end_point = '50 25 0'
    #   sort_by = x
    #   num_points = 500
    #   outputs = csv
    #   variable = 'porepressure temperature'
    # []

    # [ptsussfull_res2]
    #     type = LineValueSampler
    #     use_displaced_mesh = false
    #     start_point = '0.1 0 0'
    #     end_point = '50 0 0'
    #     sort_by = x
    #     num_points = 500
    #     outputs = csv
    #     variable = 'porepressure temperature'
    # []


    # [ptsussfull_cap1]
    #         type = LineValueSampler
    #         use_displaced_mesh = false
    #         start_point = '0.1 30 0' #5m above res
    #         end_point = '50 30 0'
    #         sort_by = x
    #         num_points = 500
    #         outputs = csv
    #         variable = 'porepressure temperature'
    # []

    # [ptsussfull_cap2]
    #     type = LineValueSampler
    #     use_displaced_mesh = false
    #     start_point = '0.1 35 0' #5m above res
    #     end_point = '50 35 0'
    #     sort_by = x
    #     num_points = 500
    #     outputs = csv
    #     variable = 'porepressure temperature'
    # []

    # [ptsussfull_cap3]
    #     type = LineValueSampler
    #     use_displaced_mesh = false
    #     start_point = '0.1 40 0' #5m above res
    #     end_point = '50 40 0'
    #     sort_by = x
    #     num_points = 500
    #     outputs = csv
    #     variable = 'porepressure temperature'
    # []
    

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
    end_time = ${time}
    dtmax = 1e6
    nl_rel_tol = 1e-4
    nl_abs_tol = 1e-8
    [TimeStepper]
      type = IterationAdaptiveDT
      dt = 1000
    []
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

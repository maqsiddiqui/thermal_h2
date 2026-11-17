########
#### THIS CODE SOLVES Lauwerier problem for H2:
########

# The mesh used below has very high resolution, so the simulation takes a long time to complete.
# Some suggested meshes of different resolution:
# nx=50, bias_x=1.2
# nx=100, bias_x=1.1
# nx=200, bias_x=1.05
# nx=400, bias_x=1.02
# nx=1000, bias_x=1.01
# nx=2000, bias_x=1.003
#####################################

inject_fluid_mass = 3.6e6 # kg # 20 ton/day equivalent


inject_time = 15552000 # # 180days
#store_time =  432000 #5*86400 # days
#produce_time = 15552000 # # days
#rest_time = 0 # days
num_cycles = 1 # Currently needs to be <= 10

cycle_length = '${fparse inject_time}'

end_simulation = '${fparse cycle_length * num_cycles}'

# Note: I have setup 10 cycles but you can set num_cycles less than 10.
start_injection1 = 0

end_injection1 = '${fparse start_injection1 + inject_time}'

plot_time1 = '${fparse end_injection1 * 0.5}'

synctimes = '${plot_time1} ${end_injection1}'

# Geometry in RZ coordinates
# borehole radius (m)
bh_r = 0.1
# model radius (m)
max_r = 1000
# aquifer thickness (m)
aq_thickness = 50
# cap thickness (m)
cap_thickness = 40
# injection region top and bottom (m).  Note, the mesh is created with the aquifer in y = (-0.5 * aq_thickness, 0.5 * aq_thickness), irrespective of depth (depth only sets the insitu porepressure and temperature)
screen_top = '${fparse 0.5 * aq_thickness}'
screen_bottom = '${fparse -0.5 * aq_thickness}'


# number of elements in radial direction
num_r = 200
# number of elements across half height of aquifer
num_y_aq = 10
# number of elements across height of cap
num_y_cap = 8

# mesh bias in radial direction
bias_r = 1.05 #1.22
# mesh bias in vertical direction in aquifer top
bias_y_aq_top = 0.9
# mesh bias in vertical direction in cap top
bias_y_cap_top = 1.3
# mesh bias in vertical direction in aquifer bottom
bias_y_aq_bottom = '${fparse 1.0 / bias_y_aq_top}'
# mesh bias in vertical direction in cap bottom
bias_y_cap_bottom = '${fparse 1.0 / bias_y_cap_top}'
#depth_centre = '${fparse depth + aq_thickness/2}'

#####################################
half_aq_thickness = '${fparse aq_thickness * 0.5}'
half_height = '${fparse half_aq_thickness + cap_thickness}'
approx_screen_length = '${fparse screen_top - screen_bottom}'

aq_density = 2601.2 # kg/m^3
aq_specific_heat_cap = 951.6 # J/Kg/K

# Thermal radius (note this is not strictly correct, it should use the bulk specific heat
#  capacity as defined below, but it doesn't matter here because this is purely for
#  defining the region of refined mesh)
th_r = '${fparse sqrt(inject_fluid_mass / 1000 * 4.12e6 / (approx_screen_length * 3.1416 * aq_specific_heat_cap * aq_density))}'
# radius of fine mesh
fine_r = '${fparse th_r * 2}'
bias_r_fine = 1
num_r_fine = '${fparse int(fine_r/1)}'



[Mesh]
    coord_type = RZ
    construct_side_list_from_node_list=true
    [aq_top_fine]
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${num_r_fine}
        xmin = ${bh_r}
        xmax = ${fine_r}
        bias_x = ${bias_r_fine}
        bias_y = ${bias_y_aq_top}
        ny = ${num_y_aq}
        ymin = 0
        ymax = ${half_aq_thickness}

    []

    [cap_top_fine]
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${num_r_fine}
        xmin = ${bh_r}
        xmax = ${fine_r}
        bias_x = ${bias_r_fine}
        bias_y = ${bias_y_cap_top}
        ny = ${num_y_cap}
        ymax = ${half_height}
        ymin = ${half_aq_thickness}
      []

    [aq_and_cap_top_fine]
        type = StitchedMeshGenerator
        inputs = 'aq_top_fine cap_top_fine'
        clear_stitched_boundary_ids = true
        stitch_boundaries_pairs = 'top bottom'
        prevent_boundary_ids_overlap = false
    []

    [aq_bottom_fine]
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${num_r_fine}
        xmin = ${bh_r}
        xmax = ${fine_r}
        bias_x = ${bias_r_fine}
        bias_y = ${bias_y_aq_bottom}
        ny = ${num_y_aq}
        ymax = 0
        ymin = -${half_aq_thickness}
    []

    [cap_bottom_fine]
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${num_r_fine}
        xmin = ${bh_r}
        xmax = ${fine_r}
        bias_x = ${bias_r_fine}
        bias_y = ${bias_y_cap_bottom}
        ny = ${num_y_cap}
        ymin = -${half_height}
        ymax = -${half_aq_thickness}
    []
      
    [aq_and_cap_bottom_fine]
        type = StitchedMeshGenerator
        inputs = 'aq_bottom_fine cap_bottom_fine'
        clear_stitched_boundary_ids = true
        stitch_boundaries_pairs = 'bottom top'
        #merge_boundaries_with_same_name = false
        prevent_boundary_ids_overlap = false
    []

    [aq_and_cap_fine]
        type = StitchedMeshGenerator
        inputs = 'aq_and_cap_bottom_fine aq_and_cap_top_fine'
        clear_stitched_boundary_ids = true
        stitch_boundaries_pairs = 'top bottom'
        prevent_boundary_ids_overlap = false
    []

    [aq_top]
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${num_r}
        xmin = ${fine_r}
        xmax = ${max_r}
        bias_x = ${bias_r}
        bias_y = ${bias_y_aq_top}
        ny = ${num_y_aq}
        ymin = 0
        ymax = ${half_aq_thickness}
    []

    [cap_top]
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${num_r}
        xmin = ${fine_r}
        xmax = ${max_r}
        bias_x = ${bias_r}
        bias_y = ${bias_y_cap_top}
        ny = ${num_y_cap}
        ymax = ${half_height}
        ymin = ${half_aq_thickness}
    []

    [aq_and_cap_top]
        type = StitchedMeshGenerator
        inputs = 'aq_top cap_top'
        clear_stitched_boundary_ids = true
        stitch_boundaries_pairs = 'top bottom'
        prevent_boundary_ids_overlap = false
    []

    [aq_bottom]
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${num_r}
        xmin = ${fine_r}
        xmax = ${max_r}
        bias_x = ${bias_r}
        bias_y = ${bias_y_aq_bottom}
        ny = ${num_y_aq}
        ymax = 0
        ymin = -${half_aq_thickness}
    []

    [cap_bottom]
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${num_r}
        xmin = ${fine_r}
        xmax = ${max_r}
        bias_x = ${bias_r}
        bias_y = ${bias_y_cap_bottom}
        ny = ${num_y_cap}
        ymin = -${half_height}
        ymax = -${half_aq_thickness}
    []

    [aq_and_cap_bottom]
        type = StitchedMeshGenerator
        inputs = 'aq_bottom cap_bottom'
        clear_stitched_boundary_ids = true
        stitch_boundaries_pairs = 'bottom top'
        prevent_boundary_ids_overlap = false
    []

    [aq_and_cap]
        type = StitchedMeshGenerator
        inputs = 'aq_and_cap_bottom aq_and_cap_top'
        clear_stitched_boundary_ids = true
        stitch_boundaries_pairs = 'top bottom'
        prevent_boundary_ids_overlap = false
    []
    
    [aq_and_cap_all]
        type = StitchedMeshGenerator
        inputs = 'aq_and_cap_fine aq_and_cap'
        clear_stitched_boundary_ids = true
        stitch_boundaries_pairs = 'right left'
        prevent_boundary_ids_overlap = false
    []
    
    [aquifer]
        type = ParsedSubdomainMeshGenerator
        input = aq_and_cap_all
        combinatorial_geometry = 'y >= -${half_aq_thickness} & y <= ${half_aq_thickness}'
        block_id = 1
    []

    [top_cap]
        type = ParsedSubdomainMeshGenerator
        input = aquifer
        combinatorial_geometry = 'y >= ${half_aq_thickness}'
        block_id = 2
    []

    [bottom_cap]
        type = ParsedSubdomainMeshGenerator
        input = top_cap
        combinatorial_geometry = 'y <= -${half_aq_thickness}'
        block_id = 3
    []

    [injection_area]
        type = ParsedGenerateSideset
        combinatorial_geometry = 'x<=${bh_r}*1.000001 & y >= ${screen_bottom} & y <= ${screen_top}'
        included_subdomains = 1
        new_sideset_name = 'injection_area'
        input = 'bottom_cap'
    []

    [rename]
        type = RenameBlockGenerator
        old_block = '1 2 3'
        new_block = 'aquifer caps caps'
        input = 'injection_area'
    []
    
[]
  

[GlobalParams]
    PorousFlowDictator = dictator
    gravity = '0 0 0'
[]

[Variables]
    [porepressure]
      initial_condition = 14.5e6 #2170 psi
    []
    [temperature]
      initial_condition = 341.15 #68 #341.15 #TVD 1470 m , 2.55C/100m
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
      value = 14.5e6
      boundary = 'bottom right top'
    #   boundary = 'right'
    []
    [outer_boundary_temperature]
      type = DirichletBC
      preset = true
      variable = temperature
    #   function = insitu_temperature
      value = 341.15
      boundary = 'bottom right top'
    #   boundary = 'right'
    []

    [inject_heat]
        type = DirichletBC
        variable = temperature
        # function = ${inject_temp}
        boundary = 'injection_area'
        value = 318.15
    []

    [inject_fluid]
      type = PorousFlowSink
      variable = porepressure
      boundary = injection_area
      flux_function = injection_rate_value
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

[]

[Functions]

    [inject]
      type = ParsedFunction
      expression = 'if(t >= ${start_injection1} & t < ${end_injection1}, 1,0)'
    
    []
  
    [injection_rate_value]
      type = ParsedFunction
      symbol_names = true_screen_area
      symbol_values = true_screen_area
      expression = '-${inject_fluid_mass}/(true_screen_area * ${inject_time})'
    []
  
[]


[FluidProperties]
    [true_water]
      type = SimpleFluidProperties
      #type = HydrogenFluidProperties #SimpleFluidProperties
      bulk_modulus = 14.53e6 #6.7e10 #2.27e14
      #density = 10 # 905.24
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
      permeability = '4.94e-13 0 0  0 4.94e-13 0  0 0 4.94e-13'
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
     
      dry_thermal_conductivity = '0 0 0  0 20 0  0 0 0'
      #wet_thermal_conductivity = '3.083 0 0  0 0 0  0 0 3.083'
      #wet_thermal_conductivity = '0 0 0  0 0 0  0 0 3.083'
      block = aquifer
    []

    [thermal_conductivity_caps]
        type = PorousFlowThermalConductivityIdeal
        block = caps
        dry_thermal_conductivity = '0 0 0  0 4.310 0  0 0 0'
        #wet_thermal_conductivity = '0 0 0  0 4.310 0  0 0 0'
    []

    [internal_energy_reservoir]
      type = PorousFlowMatrixInternalEnergy
      specific_heat_capacity = 739.95
      density = 2898.68
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
[]
  
[Kernels]
    [mass_gas_dot]
      type = PorousFlowMassTimeDerivative
      fluid_component = 0
      variable = porepressure
      multiply_by_density = true
      
    []
    [flux_gas]
      type =  PorousFlowAdvectiveFlux # PorousFlowFullySaturatedDarcyFlow
      fluid_component = 0
      use_displaced_mesh = false
      variable = porepressure
      #gravity = '0 0 0'
     
    []
    [energy_dot]
      type = PorousFlowEnergyTimeDerivative
      variable = temperature
    []
    [advection]
      type = PorousFlowHeatAdvection
      use_displaced_mesh = false
      variable = temperature
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
    [ptsussfull1]
      type = LineValueSampler
      use_displaced_mesh = false
      start_point = '0.1 10 0'
      end_point = '1000 10 0'
      sort_by = x
      num_points = 1000
      outputs = csv
      variable = 'porepressure temperature'
    []

    [ptsussfull2]
            type = LineValueSampler
            use_displaced_mesh = false
            start_point = '0.1 35 0'
            end_point = '1000 35 0'
            sort_by = x
            num_points = 1000
            outputs = csv
            variable = 'porepressure temperature'
    []

    # [ptsussfull3]
    #         type = LineValueSampler
    #         use_displaced_mesh = false
    #         start_point = '0.1 -45 0'
    #         end_point = '1000 -45 0'
    #         sort_by = x
    #         num_points = 1000
    #         outputs = csv
    #         variable = 'porepressure temperature'
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
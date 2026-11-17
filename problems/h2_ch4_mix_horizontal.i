time = 1.577e7

synctimes = '2.628e6 7.884e6 1.577e7' # 1 mo, 3mo, 6 mo

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    ny = 10
    xmin = 0
    xmax = 100
    nx= 100
    ymin = 0
    ymax = 50
  []
  [top_box]
    type = ParsedSubdomainMeshGenerator
    input = gen
    combinatorial_geometry = 'x >= 0 & x <= 15'
    block_id = 1
    block_name = top_box
  []
  [bot_box]
    type = ParsedSubdomainMeshGenerator
    input = top_box
    combinatorial_geometry = 'x >= 15 & x < 100'
    block_id = 2
    block_name = bot_box
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [pp_CO2]
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

[Kernels]
  [mass_CO2]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = pp_CO2
  []
  [adv_CO2]
    type = PorousFlowFullySaturatedDarcyFlow
    variable = pp_CO2
    fluid_component = 1
  []
  [diff_CO2]
    type = PorousFlowDispersiveFlux
    fluid_component = 1
    variable = pp_CO2
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

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp_CO2 mass_frac_CH4'
    number_fluid_phases = 1
    number_fluid_components = 2
  []
[]

[Functions]
    [s_curve]
        type = ParsedFunction
        expression = '1 / (1 + exp(-k * (x + y - m)))'
        vars = 'k m'
        vals = '5 5'
    []

[]

[FluidProperties]
  [CO2]
    type = HydrogenFluidProperties
  []
  [CH4]
    type = MethaneFluidProperties
  []
[]

[Materials]
  [ps]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp_CO2
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.18
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '6.6646e-13 0 0   0 6.6646e-13 0   0 0 6.6646e-13'
  []
  [CO2]
    type = PorousFlowMultiComponentGasMixture
    fp = 'CH4 CO2'
    mass_fraction = mass_frac_CH4
    phase = 0
    compute_enthalpy = false
    compute_internal_energy = false
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = mass_frac_CH4
  []
  [temperature]
    type = PorousFlowTemperature
    temperature = 313
  []
  [diff]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '1e-7 1e-7'
    tortuosity = 1
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 1
    phase = 0
  []
[]

[ICs]
  [mass_frac_CH4_bot]
    type = ConstantIC
    variable = mass_frac_CH4
    value = 0
    block = bot_box
  []
  [mass_frac_CH4_top]
    type = FunctionIC # ConstantIC
    variable = mass_frac_CH4
    #value = 1
    function = s_curve
    block = top_box
  []
  [pressure]
    type = ConstantIC
    variable = pp_CO2
    value = 14.5e6
  []
[]

[Preconditioning]
  [basic]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
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

[VectorPostprocessors]
  [mass_frac_CH4]
    type = LineValueSampler
    start_point = '0 50 0'
    end_point = '100 50 0'
    variable = mass_frac_CH4
    sort_by = x
    num_points = 500
  []
  [mass_frac_CH4_2]
    type = LineValueSampler
    start_point = '0 25 0'
    end_point = '100 25 0'
    variable = mass_frac_CH4
    sort_by = x
    num_points = 500
  []
  [density]
    type = LineValueSampler
    start_point = '0 50 0'
    end_point = '100 50 0'
    variable = density
    sort_by = x
    num_points = 500
  []
  [density_2]
    type = LineValueSampler
    start_point = '0 25 0'
    end_point = '100 25 0'
    variable = density
    sort_by = x
    num_points = 500
  []
[]


[Outputs]
  print_linear_residuals = false
  sync_times =  ${synctimes}
#   file_base = h2_ch4_mix_horizontal_${time}
  perf_graph = true
  [CSV]
    type = CSV
    sync_only = true
  []
  exodus = true
[]
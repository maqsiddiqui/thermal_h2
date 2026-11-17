time = 21000 # 350 mins

synctimes = '3000 6000 12000 15000' # 1 mo, 3mo, 6 mo

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    #ny = 10
    xmin = 0
    xmax = 0.20 # m
    nx= 30
    #ymin = 0
    #ymax = 50
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [pp_CO2]
  []
  [mass_frac_H2]
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

  [diff_CO2]
    type = PorousFlowDispersiveFlux
    fluid_component = 1
    variable = pp_CO2
    disp_trans = 0
    disp_long = 0
  []
  [mass_H2]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = mass_frac_H2
  []
  [diff_H2]
    type = PorousFlowDispersiveFlux
    fluid_component = 0
    variable = mass_frac_H2
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
    porous_flow_vars = 'pp_CO2 mass_frac_H2'
    number_fluid_phases = 1
    number_fluid_components = 2
  []
[]

# [Functions]
#     [s_curve]
#         type = ParsedFunction
#         expression = '1 / (1 + exp(-k * (x + y - m)))'
#         vars = 'k m'
#         vals = '5 5'
#     []

# []

[FluidProperties]
  [H2]
    type = HydrogenFluidProperties
  []
  [CO2]
    type = CO2FluidProperties
  []
[]

[Materials]
  [ps]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp_CO2
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '6.6646e-22 0 0   0 0 0   0 0 0'
  []
  [H2]
    type = PorousFlowMultiComponentGasMixture
    fp = 'CO2 H2'
    mass_fraction = mass_frac_H2
    phase = 0
    compute_enthalpy = false
    compute_internal_energy = false
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = mass_frac_H2
  []
  [temperature]
    type = PorousFlowTemperature
    temperature = 318.15
  []
  [diff]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '27.7e-8 27.7e-8'
    tortuosity = 1
  []
  [relperm]
    type = PorousFlowRelativePermeabilityConst
    # n = 1
    phase = 0
  []
[]

[ICs]
  [mass_frac_H2]
    type = ConstantIC
    variable = mass_frac_H2
    value = 0
  []

  [pressure]
    type = ConstantIC
    variable = pp_CO2
    value = 10e6
  []
[]

[BCs]

    [outer_boundary_porepressure]
        type = DirichletBC
        preset = true
        variable = pp_CO2
      #   function = insitu_pressure
        value = 10e6
        boundary = right
    []

    [inject_H2]
      type = DirichletBC
      variable = mass_frac_H2
      value = 0.5
      boundary = right
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
  dtmax = 60
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-8
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 30
  []
[]

[VectorPostprocessors]
  [mass_frac_H2]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '0.2 0 0'
    variable = mass_frac_H2
    sort_by = x
    num_points = 30
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
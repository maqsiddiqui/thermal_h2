
## S1 - base case ##

# Tmax =  353.15 # 80oC 
# Topt = 310.15 # 37 oC

inject_fluid_mass = 900000 # 900t #3614484 # kg --> 3600 t
inject_fluid_mass_h = '${fparse inject_fluid_mass*0.13/1.002}' # kg
inject_fluid_mass_c = '${fparse inject_fluid_mass*0.87/1.002}' # kg

produce_fluid_mass = '${fparse inject_fluid_mass_h + inject_fluid_mass_c}'

inject_time =  1.577e7 # # 6 months in seconds

store_time =  0 # seconds
produce_time = 1.577e7 #1.577e7 # # 6 months in seconds
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


synctimes = '${start_injection1} ${end_injection1} ${end_production1} 
${end_injection2} ${end_production2}
${end_injection3} ${end_production3}
${end_injection4} ${end_production4}
${end_injection5} ${end_production5}'


[Mesh]
    construct_side_list_from_node_list=true
    coord_type = RZ
    [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1000
    #bias_x = 1.005
    xmin = 0.1
    xmax = 1000.1
    ny = 1
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


[GlobalParams]
    PorousFlowDictator = dictator
    gravity = '0 0 0'
[]

[UserObjects]
    [dictator]
      type = PorousFlowDictator
      porous_flow_vars = 'porepressure cgH2 cgCH4 cgCO2 temperature'
      number_fluid_phases = 1
      number_fluid_components = 4
    []
    [pc]
      type = PorousFlowCapillaryPressureConst
      pc = 0
    []

    # [fluid_advective_flux]
    #     type = PorousFlowAdvectiveFluxCalculatorSaturated
    #     flux_limiter_type = superbee
    # []

    # [heat_advective_flux]
    #   type = PorousFlowAdvectiveFluxCalculatorSaturatedHeat
    #   flux_limiter_type = superbee
    # []

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

    [temperature]
      initial_condition = 336.15 #63 
      scaling = 1e-4 # 1E-6 # fluid enthalpy is roughly 1E6
    []
[]

[ICs]

    [porepressure]
        type = FunctionIC
        variable = porepressure
        function = 10e6
    []

[]
  
[AuxVariables]

     [bounds_dummy]
        order = FIRST
        family = LAGRANGE
    []

    [h2_flux_out_kg_per_s]
    []
    [co2_flux_out_kg_per_s]
    []
    [ch4_flux_out_kg_per_s]
    []

    [h2_in_kg_per_s]   
    []
    [co2_in_kg_per_s]
    []

    [h2_prod_kg_per_s]
    []
    [co2_prod_kg_per_s]
    []
    [ch4_prod_kg_per_s]
    []
    [carrier_prod_kg_per_s]
    []
[]
  
[Kernels]

    [npop_dt]
        type = TimeDerivative
        variable = npop_d
    []

    [npop_diff]
        type = MatDiffusion
        variable = npop_d
        diffusivity = 25e-6 #m2/s Dutta et al. eliFE paper

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

    [diff_gas]
        type = PorousFlowDispersiveFlux
        fluid_component = 3
        variable = porepressure
        disp_trans = 0
        disp_long = 0
    []
    # Water source from methanogenesis (4H2 + CO2 -> CH4 + 2H2O):
    # 1 + 5.5 units consumed, 2 units produced as CH4 => 4.5 units balance to water (component 3).
    # Without this term the total mass balance has a net sink, causing a spurious pressure drop
    # and an associated ~4 degC temperature artefact at the reaction zone.

     [water_from_rxn]
         type = MatCoupledForce
         variable = porepressure
         coef = 1
         v = cgH2
         material_properties = 'doublemonodterm_water'
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
    [cgCH4_diff]
        type = PorousFlowDispersiveFlux
        fluid_component = 1
        variable = cgCH4
        disp_long = 0
        disp_trans = 0
    []
    [cgCH4_conv]
        type = PorousFlowFullySaturatedDarcyFlow
        fluid_component = 1
        variable = cgCH4
        multiply_by_density = true
    []

    [cgCH4_reaction]
        type = MatCoupledForce
        variable = cgCH4
        coef = 1
        v = cgH2
        material_properties = 'doublemonodterm_ch4'
        
    []

    [cgCO2_dt]
        type = PorousFlowMassTimeDerivative
        fluid_component = 2
        variable = cgCO2
        multiply_by_density = true
    []
    [cgCO2_diff]
        type = PorousFlowDispersiveFlux
        fluid_component = 2
        variable = cgCO2
        disp_long = 0
        disp_trans = 0
    []
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


[Functions]
    [injection_rate]
      type = ParsedFunction
      symbol_values = injection_area
      symbol_names = area
      expression = '-${inject_fluid_mass}/(area* ${inject_time})' 
    []

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

    [H2_injection_rate]
      type = ParsedFunction
      symbol_values = injection_area
      symbol_names = area
    #   expression = '-${inject_fluid_mass}*0.13/(area *0.162*1.02 * ${inject_time})' #Retard_h = 1.02
      expression = '-${inject_fluid_mass_h}/(area *0.162* ${inject_time})' #Retard_h = 1.02
    []
    [CO2_injection_rate]
      type = ParsedFunction
      symbol_values = injection_area
      symbol_names = area
    #   expression = '-${inject_fluid_mass}*0.87/(area *0.162 *3.05 * ${inject_time})'
      expression = '-${inject_fluid_mass_c}/(area*0.162 * ${inject_time})'
    []

    [total_injection_rate]
      type = ParsedFunction
      symbol_values = 'H2_injection_rate CO2_injection_rate'
      symbol_names = 'h2_rate co2_rate'
      expression = 'h2_rate + co2_rate'
    []

    [production_rate]
      type = ParsedFunction
      symbol_values = injection_area
      symbol_names = area
    #   expression = '${inject_fluid_mass}*0.13/(area *0.162*1.02 * ${inject_time})' #Retard_h = 1.02
      expression = '${produce_fluid_mass}/(area *0.162 * ${produce_time})' #Retard_h = 1.02
    []

    # [h2_out_in_timestep]
    #     type = ParsedFunction
    #     symbol_names = 'dt h2_out flux_out'
    #     symbol_values = 'dt h2_kg_per_s production_flux' #these are from Postprocessor
    #     expression = 'dt*h2_out*flux_out'
    # []

    # [co2_out_in_timestep]
    #     type = ParsedFunction
    #     symbol_names = 'dt co2_out flux_out'
    #     symbol_values = 'dt co2_kg_per_s production_flux'
    #     expression = 'dt*co2_out*flux_out'
    # []

    # [ch4_out_in_timestep]
    #     type = ParsedFunction
    #     symbol_names = 'dt ch4_out flux_out'
    #     symbol_values = 'dt ch4_kg_per_s production_flux'
    #     expression = 'dt*ch4_out*flux_out'
    # []

    [h2_in_fn]
        type = ParsedFunction
        symbol_names  = 'dt r'
        symbol_values = 'dt h2_in_rate'
        expression    = 'dt*abs(r)'
    []
    [co2_in_fn]
        type = ParsedFunction
        symbol_names  = 'dt r'
        symbol_values = 'dt co2_in_rate'
        expression    = 'dt*abs(r)'
    []

    [h2_out_fn]      type = ParsedFunction symbol_names = 'dt r' symbol_values = 'dt h2_prod_rate'      expression = 'dt*abs(r)' []
    [co2_out_fn]     type = ParsedFunction symbol_names = 'dt r' symbol_values = 'dt co2_prod_rate'     expression = 'dt*abs(r)' []
    [ch4_out_fn]     type = ParsedFunction symbol_names = 'dt r' symbol_values = 'dt ch4_prod_rate'     expression = 'dt*abs(r)' []
    #[carrier_out_fn] type = ParsedFunction symbol_names = 'dt r' symbol_values = 'dt carrier_prod_rate' expression = 'dt*abs(r)' []

[]

[Controls]

    [inject_phase]
        type = ConditionalFunctionEnableControl
        enable_objects = 'BCs::H2_injection BCs::CO2_injection BCs::inject_temp'
        #enable_objects = 'BCs::total_injection BCs::cgH2_in BCs::cgCO2_in BCs::inject_temp'
        conditional_function = inject
        implicit = false
        execute_on = 'initial timestep_begin'
    []
    [produce_phase]
        type = ConditionalFunctionEnableControl
        #enable_objects = 'BCs::total_production BCs::outflow_h2 BCs::outflow_co2 BCs::outflow_ch4'
        #enable_objects = 'BCs::outflow_pp BCs::outflow_h2 BCs::outflow_co2 BCs::outflow_ch4'
        enable_objects = 'BCs::H2_production BCs::CO2_production BCs::CH4_production BCs::total_production'
        conditional_function = produce
        implicit = false
        execute_on = 'initial timestep_begin'
    []

[]

  
[FluidProperties]

    [hydrogen]
        type= IdealGasFluidProperties
    []

    # [hydrogen]
    #     type=SimpleFluidProperties
    #     bulk_modulus = 10e6 #14.5e6
    #     density0 = 9.55 # 905.24
    #     viscosity = 9.93e-6 #
    #     cv  = 1 # 4100  #10183
    #     cp = 14657.2
    #     porepressure_coefficient =0 #1e-7 #1.0
    #     thermal_expansion = 0 # 2.976e-3
    #     thermal_conductivity = 0.22
    # []

    #  [hydrogen]
    #   type = HydrogenFluidProperties # BrineFluidProperties #   HydrogenFluidProperties # BrineFluidProperties # #HydrogenFluidProperties #BrineFluidProperties
    #  []

[]
  
[Materials]


    [npop_source]
        type = ParsedMaterial
        property_name = npop_source
        constant_names = 'n0 Pgrowth_max Tmax Topt Y_h2 Sw Mg Mratio poros gamma halfsat_h2 halfsat_co2 Hc_h Hc_c Ret decay_rate'
        constant_expressions = '1.8e14 1.95e-4 353.15 310.15 4.2e11 0.1 0.002 0.022 0.18 1 2e-8 4.4e-5 0.02 1.59 1.002 9e-8'
        coupled_variables = 'cgH2 cgCO2 npop_d temperature'
        expression = 'if (cgH2>0.01, (npop_d*Pgrowth_max*(Tmax-temperature)/(Tmax-Topt)*Hc_c*Hc_h*cgCO2*cgH2/Y_h2/(halfsat_h2+(cgH2*Hc_h))/(halfsat_co2+(cgCO2*Hc_c))) - (decay_rate*npop_d),0)'
        #outputs = exodus
    []

    [doublemonodsource_h2]
        type = ParsedMaterial
        property_name = doublemonodterm_h2
        constant_names = 'n0 Pgrowth_max Tmax Topt Y_h2 Sw Mg Mratio poros gamma halfsat_h2 halfsat_co2 Hc_h Hc_c Ret'
        constant_expressions = '1.8e14 1.95e-4 353.15 310.15 4.2e11 0.1 0.002 0.022 0.18 1 2e-8 4.4e-5 0.02 1.59 1.002'
        coupled_variables = 'cgH2 cgCO2 npop_d temperature'
        expression = '-gamma*Mratio*poros*Sw*n0*npop_d*Pgrowth_max*(Tmax-temperature)/(Tmax-Topt)*Mg*Hc_c*Hc_h*cgCO2/Y_h2/Ret/(halfsat_h2+(cgH2*Hc_h))/(halfsat_co2+(cgCO2*Hc_c))'
        #outputs = exodus
    []

    [doublemonodsource_co2]
        type = ParsedMaterial
        property_name = doublemonodterm_co2
        constant_names = 'n0 Pgrowth_max Tmax Topt Y_h2 Sw Mg Mratio poros gamma halfsat_h2 halfsat_co2 Hc_h Hc_c Ret'
        constant_expressions = '1.8e14 1.95e-4 353.15 310.15 4.2e11 0.1 0.002 0.022 0.18 5.5 2e-8 4.4e-5 0.02 1.59 1.002'
        coupled_variables = 'cgH2 cgCO2 npop_d temperature'
        expression = 'gamma*Mratio*poros*Sw*n0*npop_d*Pgrowth_max*(Tmax-temperature)/(Tmax-Topt)*Mg*Hc_c*Hc_h*cgCO2/Y_h2/Ret/(halfsat_h2+(cgH2*Hc_h))/(halfsat_co2+(cgCO2*Hc_c))'
        #outputs = exodus
    []

    [doublemonodsource_ch4]
        type = ParsedMaterial
        property_name = doublemonodterm_ch4
        constant_names = 'n0 Pgrowth_max Tmax Topt Y_h2 Sw Mg Mratio poros gamma halfsat_h2 halfsat_co2 Hc_h Hc_c Ret'
        constant_expressions = '1.8e14 1.95e-4 353.15 310.15 4.2e11 0.1 0.002 0.022 0.18 2 2e-8 4.4e-5 0.02 1.59 1.002'
        coupled_variables = 'cgH2 cgCO2 npop_d temperature'
        expression = 'gamma*Mratio*poros*Sw*n0*npop_d*Pgrowth_max*(Tmax-temperature)/(Tmax-Topt)*Mg*Hc_c*Hc_h*cgCO2/Y_h2/Ret/(halfsat_h2+(cgH2*Hc_h))/(halfsat_co2+(cgCO2*Hc_c))'
        #outputs = exodus
    []

    # Mass-balance closure: water produced by methanogenesis (gamma = 1+5.5-2 = 4.5).
    # Same Monod kernel as CO2/CH4 but positive (source) and applied to the porepressure variable.
    [doublemonodsource_water]
        type = ParsedMaterial
        property_name = doublemonodterm_water
        constant_names = 'n0 Pgrowth_max Tmax Topt Y_h2 Sw Mg Mratio poros gamma halfsat_h2 halfsat_co2 Hc_h Hc_c Ret'
        constant_expressions = '1.8e14 1.95e-4 353.15 310.15 4.2e11 0.1 0.002 0.022 0.18 4.5 2e-8 4.4e-5 0.02 1.59 1.002'
        coupled_variables = 'cgH2 cgCO2 npop_d temperature'
        expression = 'gamma*Mratio*poros*Sw*n0*npop_d*Pgrowth_max*(Tmax-temperature)/(Tmax-Topt)*Mg*Hc_c*Hc_h*cgCO2/Y_h2/Ret/(halfsat_h2+(cgH2*Hc_h))/(halfsat_co2+(cgCO2*Hc_c))'
    []


    [temperature]
      type = PorousFlowTemperature
      temperature = temperature #
      #at_nodes = true
    []
    [ppss]
        type = PorousFlow1PhaseFullySaturated
        porepressure = porepressure
    []

    [massfrac]
        type = PorousFlowMassFraction
        mass_fraction_vars = 'cgH2 cgCH4 cgCO2'
    []

    [hydrogen]
      type =  PorousFlowSingleComponentFluid # PorousFlowMultiComponentFluid # PorousFlowSingleComponentFluid #
      fp = hydrogen
      #fp = H2
      #x = xH2
      phase = 0
      #at_nodes = true
    []

    [porosity_reservoir]
      type = PorousFlowPorosityConst
      porosity = 0.162
      #at_nodes = true
    []

    [permeability_reservoir]
      type = PorousFlowPermeabilityConst
      #permeability = '1.8e-11 0 0  0 1.8e-11 0  0 0 1.8e-11' #3e-13 = 300 mD
      permeability = '6.6646e-13 0 0  0 6.6646e-13 0  0 0 0'
    []

    [relperm]
     type = PorousFlowRelativePermeabilityConst
     phase = 0
     #at_nodes = true
     s_res= 0.1
     sum_s_res = 0.1
    []

    [diff]
        type = PorousFlowDiffusivityConst
        diffusion_coeff = '8.94e-7 8.94e-7 8.94e-7 8.94e-7'
        tortuosity = 1
    []

    [thermal_conductivity_reservoir]
      type = PorousFlowThermalConductivityIdeal
     
      dry_thermal_conductivity = '3.083 0 0  0 20 0 0 0 0'
      #wet_thermal_conductivity = '3.083 0 0  0 0 0  0 0 3.083'
      #wet_thermal_conductivity = '0 0 0  0 0 0  0 0 3.083'
    []

    [internal_energy_reservoir]
      type = PorousFlowMatrixInternalEnergy
      specific_heat_capacity = 739.95
      density = 2898.68
    []

  []
  
[BCs]


    [H2_injection]
      type = PorousFlowSink
      boundary = injection_area
      variable = cgH2
      use_mobility =  false
      use_relperm = false
      flux_function =  H2_injection_rate
      #value = -0.000125
    #   mass_fraction_component = 0
    #   fluid_phase=0
        save_in = h2_in_kg_per_s
    []

    [CO2_injection]
      type = PorousFlowSink
      boundary = injection_area
      variable = cgCO2
      use_mobility =false
      use_relperm = false
      flux_function =  CO2_injection_rate
      #value = -0.000125
    #   mass_fraction_component = 2
    #   fluid_phase=0
        save_in = co2_in_kg_per_s
    []

    [total_production]
        type = PorousFlowPiecewiseLinearSink
        boundary = injection_area
        variable = porepressure
        fluid_phase = 0
        mass_fraction_component = 3
        use_mobility = true
        use_relperm = true
        PT_shift = 9.98e6
        pt_vals     = '0 1e7'
        multipliers = '0 1e7'
        save_in = carrier_prod_kg_per_s
    []

    [H2_production]
        type = PorousFlowPiecewiseLinearSink
        boundary = injection_area
        variable = cgH2
        fluid_phase = 0
        mass_fraction_component = 0
        use_mobility = true
        use_relperm = true
        PT_shift = 9.98e6
        pt_vals     = '0 1e7'
        multipliers = '0 1e7'
        save_in = h2_prod_kg_per_s
    []

    [CO2_production]
        type = PorousFlowPiecewiseLinearSink
        boundary = injection_area
        variable = cgCO2
        fluid_phase = 0
        mass_fraction_component = 3
        use_mobility = true
        use_relperm = true
        PT_shift = 9.98e6
        pt_vals     = '0 1e7'
        multipliers = '0 1e7'
        save_in = co2_prod_kg_per_s
    []

    [CH4_production]
        type = PorousFlowPiecewiseLinearSink
        boundary = injection_area
        variable = cgCH4
        fluid_phase = 0
        mass_fraction_component = 1
        use_mobility = true
        use_relperm = true
        PT_shift = 9.98e6
        pt_vals     = '0 1e7'
        multipliers = '0 1e7'
        save_in = ch4_prod_kg_per_s
    []


    # [outflow_pp]
    #     type = PorousFlowOutflowBC
    #     boundary = injection_area
    #     include_relperm = true #false # no need for relperm in this fully-saturated simulation
    #     mass_fraction_component = 3
    #     variable = porepressure
        
    # []

    # [outflow_h2]
    #     type = PorousFlowOutflowBC
    #     boundary = injection_area
    #     include_relperm = true # no need for relperm in this fully-saturated simulation
    #     mass_fraction_component = 0
    #     variable = cgH2
    #     save_in = h2_flux_out_kg_per_s
        
    # []

    # [outflow_co2]
    #     type = PorousFlowOutflowBC
    #     boundary = injection_area
    #     include_relperm = true # no need for relperm in this fully-saturated simulation
    #     mass_fraction_component = 2
    #     variable = cgCO2
    #     save_in = co2_flux_out_kg_per_s
        
    # []

    # [outflow_ch4]
    #     type = PorousFlowOutflowBC
    #     boundary = injection_area
    #     include_relperm = true # no need for relperm in this fully-saturated simulation
    #     mass_fraction_component = 1
    #     variable = cgCH4
    #     save_in = ch4_prod_kg_per_s 
        
    # []


    [outer_boundary_temperature]
      type = DirichletBC
      preset = true
      variable = temperature
    #   function = insitu_temperature
      value = 336.15
      boundary = 'right'
    []

    # [outer_boundary_pressure]
    #     type = FunctionDirichletBC
    #     preset = true
    #     variable = porepressure
    #     function = '10e6-x'
    #     boundary = 'right'
    # []

    [outer_boundary_pressure]
        type = DirichletBC
        preset = true
        variable = porepressure
        value = 10e6
        boundary = 'right'
    []

    [inject_temp]
        type = DirichletBC
        variable = temperature
        preset = true
        # function = ${inject_temp}
        boundary = 'injection_area'
        value = 318.15
    []

    # [produce_heat]
    #     type = PorousFlowSink
    #     variable = temperature
    #     boundary = injection_area
    #     flux_function = production_rate
    #     fluid_phase = 0
    #     use_enthalpy = true
    #     #save_in = heat_flux_out
    #     #mass_fraction_component = 3
    # []

    [produce_heat]
        type = PorousFlowOutflowBC
        variable = temperature
        boundary = injection_area
        flux_type = heat
    []


[]

[Bounds]
  [h2_lower]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = cgH2
    bound_type = lower
    bound_value = 0 #12
  []

  [co2_lower]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = cgCO2
    bound_type = lower
    bound_value =0 
  []

  [ch4_lower]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = cgCH4
    bound_type = lower
    bound_value =0 #e-12
  []


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

    [production_flux]
        type = FunctionValuePostprocessor
        function = production_rate
        execute_on = initial
    []


    [dt]
        type = TimestepSize
    []

    # [h2_kg_per_s]
    #     type = NodalSum
    #     boundary = injection_area
    #     variable = h2_flux_out_kg_per_s
    #     execute_on = 'initial timestep_end'
    #     outputs = 'none'
    # []

    # [co2_kg_per_s]
    #     type = NodalSum
    #     boundary = injection_area
    #     variable = co2_flux_out_kg_per_s
    #     execute_on = 'initial timestep_end'
    #     outputs = 'none'
    # []
    # [ch4_kg_per_s]
    #     type = NodalSum
    #     boundary = injection_area
    #     variable = ch4_flux_out_kg_per_s
    #     execute_on = 'initial timestep_end'
    #     outputs = 'none'
    # []

    # [h2_out_per_timestep]
    #     type = FunctionValuePostprocessor
    #     function = h2_out_in_timestep
    #     execute_on = 'timestep_end'
    #     outputs = 'none'
    # []

    # [h2_out_cumulative]
    #     type = CumulativeValuePostprocessor
    #     postprocessor = h2_out_per_timestep
    #     execute_on = 'timestep_end'
    #     outputs = 'csv console'
    # []

    # [co2_out_per_timestep]
    #     type = FunctionValuePostprocessor
    #     function = co2_out_in_timestep
    #     execute_on = 'timestep_end'
    #     outputs = 'none'
    # []

    # [co2_out_cumulative]
    #     type = CumulativeValuePostprocessor
    #     postprocessor = co2_out_per_timestep
    #     execute_on = 'timestep_end'
    #     outputs = 'csv console'
    # []

    # [ch4_out_per_timestep]
    #     type = FunctionValuePostprocessor
    #     function = ch4_out_in_timestep
    #     execute_on = 'timestep_end'
    #     outputs = 'none'
    # []

    # [ch4_out_cumulative]
    #     type = CumulativeValuePostprocessor
    #     postprocessor = ch4_out_per_timestep
    #     execute_on = 'timestep_end'
    #     outputs = 'csv console'
    # []

    [h2_fluid_mass]
        type = PorousFlowFluidMass
        fluid_component = 0
        execute_on = 'initial timestep_end'
        outputs = 'csv'
    []

    [ch4_fluid_mass]
        type = PorousFlowFluidMass
        fluid_component = 1
        execute_on = 'initial timestep_end'
        outputs = 'csv'
    []

    [co2_fluid_mass]
        type = PorousFlowFluidMass
        fluid_component = 2
        execute_on = 'initial timestep_end'
        outputs = 'csv'
    []

    [carrier_fluid_mass]
        type = PorousFlowFluidMass
        fluid_component = 3
        execute_on = 'initial timestep_end'
        outputs = 'csv'
    []

    [total_fluid_mass]
        type = PorousFlowFluidMass
        execute_on = 'initial timestep_end'
        outputs = 'csv'
    []

    [total_fluid_mass_initial]
        type = PorousFlowFluidMass
        execute_on = 'initial'
        outputs = 'csv'

    []

    [h2_in_rate]
        type = NodalSum
        boundary = injection_area
        variable = h2_in_kg_per_s
        execute_on = 'initial timestep_end'
        outputs = 'none'
    []
    [co2_in_rate]
        type = NodalSum
        boundary = injection_area
        variable = co2_in_kg_per_s
        execute_on = 'initial timestep_end'
        outputs = 'none'
    []
    [h2_in_cumulative]
        type = CumulativeValuePostprocessor
        postprocessor = h2_in_step
    []
    [co2_in_cumulative]
        type = CumulativeValuePostprocessor
        postprocessor = co2_in_step
    []
    [h2_in_step]
        type = FunctionValuePostprocessor
        function = h2_in_fn
        execute_on = 'timestep_end'
        outputs = 'none'
    []
    [co2_in_step]
        type = FunctionValuePostprocessor
        function = co2_in_fn
        execute_on = 'timestep_end'
        outputs = 'none'
    []

    [h2_prod_rate]      type = NodalSum boundary = injection_area variable = h2_prod_kg_per_s      execute_on = 'initial timestep_end' outputs = 'none' []
    [co2_prod_rate]     type = NodalSum boundary = injection_area variable = co2_prod_kg_per_s     execute_on = 'initial timestep_end' outputs = 'none' []
    [ch4_prod_rate]     type = NodalSum boundary = injection_area variable = ch4_prod_kg_per_s     execute_on = 'initial timestep_end' outputs = 'none' []
    #[carrier_prod_rate] type = NodalSum boundary = injection_area variable = carrier_prod_kg_per_s execute_on = 'initial timestep_end' outputs = 'none' []

    # step masses (kg) = dt * rate
    [h2_out_step]      type = FunctionValuePostprocessor function = h2_out_fn      execute_on = 'timestep_end' outputs = 'none' []
    [co2_out_step]     type = FunctionValuePostprocessor function = co2_out_fn     execute_on = 'timestep_end' outputs = 'none' []
    [ch4_out_step]     type = FunctionValuePostprocessor function = ch4_out_fn     execute_on = 'timestep_end' outputs = 'none' []
    #[carrier_out_step] type = FunctionValuePostprocessor function = carrier_out_fn execute_on = 'timestep_end' outputs = 'none' []

    # cumulative produced (kg). These REPLACE your old h2/co2/ch4_out_cumulative,
    # which were driven by the (near-zero) OutflowBCs.
    [h2_out_cumulative]
        type = CumulativeValuePostprocessor
        postprocessor = h2_out_step 
        outputs = 'csv'     
    []
    [co2_out_cumulative]     type = CumulativeValuePostprocessor postprocessor = co2_out_step     outputs = 'csv' []
    [ch4_out_cumulative]     type = CumulativeValuePostprocessor postprocessor = ch4_out_step     outputs = 'csv' []
    #[carrier_out_cumulative] type = CumulativeValuePostprocessor postprocessor = carrier_out_step outputs = 'csv' []

    # [run_time]
    #     type = PerfGraphData
    #     section_name = Root
    #     data_type = TOTAL
    # []
    

[]
  
[VectorPostprocessors]
    [ptsuss_ig]
      type = LineValueSampler
      use_displaced_mesh = false
      start_point = '0.1 2.5 0'
      end_point = '1000.1 2.5 0'
      sort_by = x
      num_points = 1000
      outputs = csv
      variable = 'cgH2 cgCH4 cgCO2 npop_d temperature'
    []

    [point_sample_ig]
          type = PointValueSampler
          variable = 'cgH2 cgCH4 cgCO2 npop_d temperature'
          points = '0.1 0 0'
          sort_by = x
          contains_complete_history = true
    []
[]
  
[Preconditioning]
    active = 'mumps'
    [smp]
      type = SMP
      full = true
      #petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
      petsc_options_iname = '-snes_type -ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
      petsc_options_value = 'vinewtonrsls gmres      asm      lu           NONZERO                   2               1E-5       1E-5        500'
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
    dtmax = 1e4 #2.5e3
    line_search = none
    [TimeStepper]
        type = IterationAdaptiveDT
        dt = 1000
        growth_factor = 2
    []
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
package Assignment_J57
  model Design
    inner BasicAeroEngines.Components.Environment environment(onDesignInit = true, Pb(displayUnit = "kPa") = 82705, Tb(displayUnit = "K") = 253.03, Mach = 0, useMach = true) annotation(
      Placement(transformation(origin = {50, 50}, extent = {{-10, -10}, {10, 10}})));
    //Intake
    BasicAeroEngines.Components.AirIntake airIntake annotation(
      Placement(transformation(origin = {-90, -38}, extent = {{-16, -16}, {16, 16}})));
    //Compressors
    BasicAeroEngines.Components.CompressorMapsBetaLines LPC(data = LPC_map, P_E(displayUnit = "kPa"), P_L(displayUnit = "kPa")) annotation(
      Placement(transformation(origin = {-60, -38}, extent = {{-10, -10}, {10, 10}})));
    BasicAeroEngines.Components.CompressorBleed HPC(data = HPC_map, Nbleed = 2, Nstages = 6, Nstages_Bleeds = {3, 4}) annotation(
      Placement(transformation(origin = {-32, -14}, extent = {{-10, -10}, {10, 10}})));
    //Combustor
    BasicAeroEngines.Components.CombustionChamberLHV combustionChamberLHV(LHV(displayUnit = "MJ/kg") = 4.296e7, V = 0.1, P_start(displayUnit = "kPa") = 1.120528e6, T_start(displayUnit = "K") = 598.4, ZC = 1, ZH = 1.9167, eta_comb = 0.9945, P_Loss = 5.5, steadyStateInit = true) annotation(
      Placement(transformation(origin = {0, 14}, extent = {{-10, -10}, {10, 10}})));
    //Turbines
    BasicAeroEngines.Components.CooledTurbine cooledHPT(data = HPT_map, eta_mech = 0.982, Nstages = 1, Nbleed = 1, Nstages_Bleeds = {1}, Xi_cool_stat = {0.05}, Xi_cool_rot = {0.05}, CoolingTechStat = {2}, CoolingTechRot = {2}) annotation(
      Placement(transformation(origin = {32, -14}, extent = {{-10, -10}, {10, 10}})));
    BasicAeroEngines.Components.CooledTurbine LPT(data = LPT_map, eta_mech = 1, Nstages = 2, CoolingTechStat = {2, 2}, CoolingTechRot = {2, 2}, Nbleed = 1, Nstages_Bleeds = {2}, Xi_cool_stat = {0, 0.02}) annotation(
      Placement(transformation(origin = {50, -38}, extent = {{-10, -10}, {10, 10}})));
  //Exhaust
    //Shafts
    BasicAeroEngines.Components.ShaftInertia HP_shaft(J = 74, omega_nom(displayUnit = "rpm") = 967.0869385300581) annotation(
      Placement(transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}})));
    BasicAeroEngines.Components.ShaftInertia LP_shaft(J = 380, omega_nom(displayUnit = "rpm") = 623.606141737574) annotation(
      Placement(transformation(origin = {0, -38}, extent = {{-10, -10}, {10, 10}})));
    //Maps
    parameter BasicAeroEngines.Data.Compressors.GSP_LPC LPC_map(P_L_nom(displayUnit = "kPa") = 326684, T_E_nom(displayUnit = "K"), P_E_nom(displayUnit = "kPa"), omega_nom(displayUnit = "rpm"), eta_nom = 0.82/0.96250, beta_nom = 11) annotation(
      Placement(transformation(origin = {-58, -64}, extent = {{-10, -10}, {10, 10}})));
    parameter BasicAeroEngines.Data.Compressors.GSP_HPC HPC_map(P_L_nom(displayUnit = "kPa") = 1.120528e6, P_E_nom(displayUnit = "kPa"), T_E_nom(displayUnit = "K"), omega_nom(displayUnit = "rpm"), eta_nom = 0.85) annotation(
      Placement(transformation(origin = {-40, 12}, extent = {{-10, -10}, {10, 10}})));
    parameter BasicAeroEngines.Data.Turbines.GSP_HPT_BetaLines HPT_map(P_L_nom(displayUnit = "kPa") = 462389, P_E_nom(displayUnit = "kPa"), T_E_nom(displayUnit = "K"), omega_nom(displayUnit = "rpm"), eta_nom = 0.865) annotation(
      Placement(transformation(origin = {44, 12}, extent = {{-10, -10}, {10, 10}})));
    parameter BasicAeroEngines.Data.Turbines.GSP_LPT_BetaLines LPT_map(P_L_nom(displayUnit = "kPa") = 237732, P_E_nom(displayUnit = "kPa"), T_E_nom(displayUnit = "K"), omega_nom(displayUnit = "rpm"), eta_nom = 0.882) annotation(
      Placement(transformation(origin = {60, -70}, extent = {{-10, -10}, {10, 10}})));
    //Fuel input
    Modelica.Blocks.Sources.TimeTable fuelFlow(table = [0, 1.03684; 10, 1.03684]) annotation(
      Placement(transformation(origin = {-16, 42}, extent = {{-10, -10}, {10, 10}})));
    //Expressions
    Modelica.Blocks.Sources.RealExpression netThrust(y = nozzleExhaustChokeable.thrust - airIntake.drag) annotation(
      Placement(transformation(origin = {142, -18}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.RealExpression thermalEfficiency(y = (nozzleExhaustChokeable.W - airIntake.W)/combustionChamberLHV.Q) annotation(
      Placement(transformation(origin = {142, -38}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.RealExpression overalEfficiency(y = (netThrust.y*environment.v)/combustionChamberLHV.Q) annotation(
      Placement(transformation(origin = {142, -60}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.RealExpression SFC(y = fuelFlow.y/nozzleExhaustChokeable.thrust) annotation(
      Placement(transformation(origin = {142, 0}, extent = {{-10, -10}, {10, 10}})));
    // Ducts
    BasicAeroEngines.Components.LinearPressureDropAir linearPressureDropAir(referenceMassFlowRate = 70.955, referencePressureDrop(displayUnit = "kPa") = 6533) annotation(
      Placement(transformation(origin = {-54, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    BasicAeroEngines.Components.BleedAirDistributor bleedAirDistributor(Nbleed = 2, useOverBleed = true, useHandBleed = false, NPortsOverBleed = {1}, NPortsLPTBleed = {2}, OverBleedPorts = 1, LPTBleedPorts = 1) annotation(
      Placement(transformation(origin = {4, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    BasicAeroEngines.Components.FlowSourceAir flowSourceAir(referenceMassFlowRate = -0.005*70.955, referenceTemperature(displayUnit = "K") = 401.05) annotation(
      Placement(transformation(origin = {18, -62}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
    BasicAeroEngines.Components.LinearPressureDropExhaust linearPressureDropExhaust(referenceMassFlowRate = 71.636, referencePressureDrop(displayUnit = "kPa") = 10572) annotation(
      Placement(transformation(origin = {68, -28}, extent = {{-6, -6}, {6, 6}})));
  Components.NozzleExhaustChokeable nozzleExhaustChokeable(f_nom = 71.636, A_fixed = 0.2375)  annotation(
      Placement(transformation(origin = {89, -37}, extent = {{-15, -15}, {15, 15}})));
  equation
    connect(airIntake.outlet, LPC.inlet) annotation(
      Line(points = {{-80, -28}, {-66, -28}}, color = {170, 223, 255}));
    connect(LP_shaft.flange_a, LPC.shaft_b) annotation(
      Line(points = {{-10, -38}, {-54, -38}}));
    connect(cooledHPT.shaft, HP_shaft.flange_b) annotation(
      Line(points = {{26, -14}, {10, -14}}));
    connect(cooledHPT.inlet, combustionChamberLHV.exhaust) annotation(
      Line(points = {{26, -8}, {26, 14}, {10, 14}}, color = {129, 170, 194}));
    connect(fuelFlow.y, combustionChamberLHV.fuelFlow) annotation(
      Line(points = {{-4, 42}, {0, 42}, {0, 24}}, color = {0, 0, 127}));
//annotation(
//  __OpenModelica_commandLineOptions = "--tearingStrictness=casual");
    connect(LPC.outlet, linearPressureDropAir.inlet) annotation(
      Line(points = {{-54, -32}, {-54, -26}}, color = {170, 223, 255}));
    connect(HPC.shaft_b, HP_shaft.flange_a) annotation(
      Line(points = {{-26, -14}, {-10, -14}}));
    connect(HPC.outlet, combustionChamberLHV.airInlet) annotation(
      Line(points = {{-26, -8}, {-26, 14}, {-10, 14}}, color = {170, 223, 255}));
    connect(HPC.inlet, linearPressureDropAir.outlet) annotation(
      Line(points = {{-38, -4}, {-54, -4}, {-54, -6}}, color = {170, 223, 255}));
    connect(HPC.outlet, cooledHPT.Bl_port[1]) annotation(
      Line(points = {{-26, -8}, {-26, -26}, {32, -26}, {32, -6}}, color = {170, 223, 255}));
    connect(LPT.shaft, LP_shaft.flange_b) annotation(
      Line(points = {{44, -38}, {10, -38}}));
    connect(LPT.inlet, cooledHPT.outlet) annotation(
      Line(points = {{44, -32}, {44, -4}, {38, -4}}, color = {129, 170, 194}));
    connect(HPC.Bl_port, bleedAirDistributor.Bl_port) annotation(
      Line(points = {{-32, -6}, {-32, -60}, {-4, -60}}, color = {170, 223, 255}, thickness = 0.5));
    connect(bleedAirDistributor.LPTBleed[1], LPT.Bl_port[1]) annotation(
      Line(points = {{4, -58}, {50, -58}, {50, -30}}, color = {170, 223, 255}, thickness = 0.5));
    connect(LPT.outlet, linearPressureDropExhaust.inlet) annotation(
      Line(points = {{56, -28}, {62, -28}}, color = {129, 170, 194}));
    connect(bleedAirDistributor.OverBleed[1], flowSourceAir.fluidPort) annotation(
      Line(points = {{4, -62}, {12, -62}}, color = {170, 223, 255}));
  connect(linearPressureDropExhaust.outlet, nozzleExhaustChokeable.inlet) annotation(
      Line(points = {{74, -28}, {80, -28}, {80, -29}}, color = {129, 170, 194}));
    annotation(
      Icon(graphics = {Text(origin = {1, 1}, extent = {{-81, 79}, {81, -79}}, textString = "ADP")}),
      Diagram,
      experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.02));
  end Design;

  model OffDesign
    extends Assignment_J57.Design(environment(onDesignInit = false, Mach = 0, useMach = true), LPC_map(P_E_nom = 83963, T_E_nom = 254.12, f_nom = 70.954, omega_nom = 623.606141737574), linearPressureDropAir(referenceMassFlowRate = 70.954), HPC_map(P_E_nom = 320151, T_E_nom = 400.35, f_nom = 70.954, omega_nom = 967.0869385300581), combustionChamberLHV(T_start = 597.11, steadyStateInit = false), HPT_map(P_E_nom = 1.058899e6, P_L_nom = 458175, T_E_nom = 1189.39, f_nom = 63.847, omega_nom = 967.0869385300581), LPT_map(P_E_nom = 458175, P_L_nom = 237108, T_E_nom = 962.89, f_nom = 70.23, omega_nom = 623.606141737574), flowSourceAir(referenceMassFlowRate = -0.354775, referenceTemperature = 400.35), fuelFlow(table = [0, 1.03684; 100, 1.03684]), nozzleExhaustChokeable(A_fixed = 0.23024082596048073));
  equation

    annotation(
      Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-100, 100}, {100, -100}}, textString = "OFF
DES")}),
  experiment(StartTime = 0, StopTime = 250, Tolerance = 1e-06, Interval = 0.5));
  end OffDesign;

  package Components
    model NozzleExhaustChokeable "Nozzle exhaust model converting static pressure into thrust"
      outer BasicAeroEngines.Components.Environment environment;
      package ExhaustFlow = BasicAeroEngines.Media.ExhaustGas;
      Modelica.Blocks.Interfaces.RealOutput thrust(unit = "N") annotation(
        Placement(visible = true, transformation(origin = {62, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, 3.55271e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      final parameter Modelica.Units.SI.Area A(fixed = false) "Nozzle outlet area";
      //parameter Modelica.SIunits.PerUnit PRnom = 1.5;
      parameter Modelica.Units.SI.Velocity v_start = 250 "Tentative value of exhaust gas velocity";
      //parameter Boolean OnDesign = false;
      parameter Modelica.Units.SI.MassFlowRate f_nom "Nominal mass flow rate" annotation(
        Dialog(enable = environment.onDesignInit));
      parameter Modelica.Units.SI.Area A_fixed "Outlet area of nozzle if known a priori";
      Modelica.Units.SI.MassFlowRate f_E "Entering mass flow rate";
      Modelica.Units.SI.Pressure P_E "Entering pressure";
      Modelica.Units.SI.Pressure P_L "Leaving pressure";
      Modelica.Units.SI.PerUnit PR_crit(start = 2.1, fixed = false) "Pressure in sonic throat";
      Modelica.Units.SI.SpecificEnthalpy h_E "Entering specific enthalpy";
      Modelica.Units.SI.SpecificEnthalpy X_E[ExhaustFlow.nX] "Entering gas composition";
      Modelica.Units.SI.SpecificEnthalpy X_L[ExhaustFlow.nX] "Entering gas composition";
      Modelica.Units.SI.SpecificEntropy s "Specific entropy of exhaust gas";
      Modelica.Units.SI.Velocity v(start = v_start) "Velocity at the nozzle exhaust";
      Modelica.Units.SI.Velocity c "Velocity of sound at the nozzle exhaust conditions";
      Modelica.Units.SI.Power W "Mechanical power output";
      ExhaustFlow.BaseProperties inletProps(T(start = 338.15)) "Fluid properties at the nozzle inlet";
      ExhaustFlow.BaseProperties outletProps(T(start = 288.15)) "Fluid properties at the nozzle outlet";
      ExhaustFlow.BaseProperties outletPropsTotal(T(start = 288.15)) "Fluid properties at the nozzle outlet";
      
      Modelica.Units.SI.PerUnit gamma_mean;
      BasicAeroEngines.Interfaces.ExhaustPort inlet annotation(
        Placement(transformation(extent = {{-68, 46}, {-48, 66}})));
    initial equation
      if environment.onDesignInit then
        f_E = f_nom;
      else
        A = A_fixed;
      end if;
    equation
// Nozzle inlet conditions (stagnation assumed)
      inletProps.p = P_E;
      inletProps.h = h_E;
      inletProps.X = X_E;
      s = ExhaustFlow.specificEntropy(inletProps.state);
// Nozzle outlet conditions
      outletProps.X = X_E;
      inletProps.h = outletProps.h + v^2/2;
      s = ExhaustFlow.specificEntropy(outletProps.state);
      c = ExhaustFlow.velocityOfSound(outletProps.state);
      gamma_mean = 0.5*(ExhaustFlow.specificHeatCapacityCp(inletProps.state)/ExhaustFlow.specificHeatCapacityCv(inletProps.state) + ExhaustFlow.specificHeatCapacityCp(outletProps.state)/ExhaustFlow.specificHeatCapacityCv(outletProps.state));
      PR_crit = (2/(gamma_mean + 1))^(-gamma_mean/(gamma_mean - 1));
      outletProps.d*v*A = f_E;
      outletProps.p = P_L;
// Nozzle outlet conditions
      if (P_E/environment.P < PR_crit) then
        P_L = environment.P;
      else
        P_L = P_E/PR_crit;
      end if;
      assert(v < c, "Invalid supersonic conditions at nozzle outlet (not a Laval nozzle)", AssertionLevel.warning);
// Generated thrust and power
      thrust = f_E*v + (outletProps.p - environment.P)*A;
      W = f_E*v^2/2;
// Boundary conditions
      P_E = inlet.P;
      f_E = inlet.f;
      h_E = inStream(inlet.h_L);
      X_E = inStream(inlet.X_L);
//X_E = {0.768,0.232};
      X_L = X_E;
      inlet.h_L = 0 "Not used, no flow reversal";
      inlet.X_L = ExhaustFlow.reference_X "Not used, no flow reversal";
// Calculating the total properties
      outletPropsTotal.h = inletProps.h;
      outletPropsTotal.X = X_E;
      outletPropsTotal.p = P_E;
      
      annotation(
        Icon(graphics = {Polygon(fillColor = {150, 150, 150}, fillPattern = FillPattern.Solid, points = {{-60, 60}, {-60, -60}, {60, -100}, {60, 100}, {60, 100}, {-60, 60}})}, coordinateSystem(initialScale = 0.1)));
    end NozzleExhaustChokeable;
  end Components;
  
  model OffDesignProperlyScaledMaps
    extends Design(environment(onDesignInit = false, Mach = 0, useMach = true), LPC_map(P_E_nom = 82705, T_E_nom = 253.03, f_nom = 70.954, omega_nom = 623.606141737574, P_L_nom = 326684, eta_nom = 0.82), linearPressureDropAir(referenceMassFlowRate = 70.954), HPC_map(P_E_nom = 320151, T_E_nom = 400.59, f_nom = 70.954, omega_nom = 967.0869385300581, P_L_nom(displayUnit = "kPa") = 1.12053e6, eta_nom = 0.85), combustionChamberLHV(T_start = 597.11, steadyStateInit = false), HPT_map(P_E_nom = 1.0589e6, P_L_nom (displayUnit = "kPa")= 458033, T_E_nom = 1189.68, f_nom = 63.85, omega_nom = 967.0869385300581, eta_nom = 0.865), LPT_map(P_E_nom = 458033, P_L_nom = 235536, T_E_nom = 962.9, f_nom = 70.23, omega_nom = 623.606141737574, eta_nom = 0.882), flowSourceAir(referenceMassFlowRate = -0.354775, referenceTemperature = 400.35), fuelFlow(table = [0, 1.03684; 100, 1.03684]), nozzleExhaustChokeable(A_fixed = 0.23024082596048073));
  equation
  
    annotation(
      Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-100, 100}, {100, -100}}, textString = "OFF
  DES")}),
      experiment(StartTime = 0, StopTime = 250, Tolerance = 1e-06, Interval = 0.5));
  end OffDesignProperlyScaledMaps;
  
  model OffDesign_Case1
    extends Design(environment(onDesignInit = false, Mach = 0, useMach = true), LPC_map(P_E_nom = 82705, T_E_nom = 253.03, f_nom = 70.954, omega_nom = 623.606141737574, P_L_nom = 326684, eta_nom = 0.82), linearPressureDropAir(referenceMassFlowRate = 70.954), HPC_map(P_E_nom = 320151, T_E_nom = 400.59, f_nom = 70.954, omega_nom = 967.0869385300581, P_L_nom(displayUnit = "kPa") = 1.12053e6, eta_nom = 0.85), combustionChamberLHV(T_start = 597.11, steadyStateInit = false), HPT_map(P_E_nom = 1.0589e6, P_L_nom (displayUnit = "kPa")= 458033, T_E_nom = 1189.68, f_nom = 63.85, omega_nom = 967.0869385300581, eta_nom = 0.865), LPT_map(P_E_nom = 458033, P_L_nom = 235536, T_E_nom = 962.9, f_nom = 70.23, omega_nom = 623.606141737574, eta_nom = 0.882), flowSourceAir(referenceMassFlowRate = -0.354775, referenceTemperature = 400.35), fuelFlow(table = [0, 0.535; 50, 0.535; 150, 1.03684]), nozzleExhaustChokeable(A_fixed = 0.23024082596048073), HP_shaft(J = 1), LP_shaft(J = 1));
  equation
  
    annotation(
      Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-100, 100}, {100, -100}}, textString = "OFF
  DES")}),
      experiment(StartTime = 0, StopTime = 250, Tolerance = 1e-06, Interval = 0.5));
  end OffDesign_Case1;

  model OnDesign_new
    extends Assignment_J57.Design(LPC(N_n_design = 1.2), LPC_map(eta_nom = 0.82/0.88750, P_L_nom = 326684, beta_nom = 15), nozzleExhaustChokeable(f_nom = 71.636, A_fixed));
  equation

  annotation(
      experiment(StartTime = 0, StopTime = 250, Tolerance = 1e-06, Interval = 0.5));
end OnDesign_new;
  
  model OffDesign_Case1_new
    extends Design(environment(onDesignInit = false, Mach = 0, useMach = true), LPC_map(P_E_nom = 82705, T_E_nom = 253.03, f_nom = 77.02684821397382, omega_nom (displayUnit = "rad/s")= 519.6717847813117, P_L_nom = 326684, eta_nom = 0.923943661971831), linearPressureDropAir(referenceMassFlowRate = 70.954), HPC_map(P_E_nom = 320151, T_E_nom = 400.59, f_nom = 70.954, omega_nom = 967.0869385300581, P_L_nom(displayUnit = "kPa") = 1.12053e6, eta_nom = 0.85), combustionChamberLHV(T_start = 597.11, steadyStateInit = false), HPT_map(P_E_nom = 1.0589e6, P_L_nom(displayUnit = "kPa") = 458033, T_E_nom = 1189.68, f_nom = 63.85, omega_nom = 967.0869385300581, eta_nom = 0.865), LPT_map(P_E_nom = 458033, P_L_nom = 235536, T_E_nom = 962.9, f_nom = 70.23, omega_nom = 623.606141737574, eta_nom = 0.882), flowSourceAir(referenceMassFlowRate = -0.354775, referenceTemperature = 400.35), fuelFlow(table = [0, 0.535; 50, 0.535; 150, 1.03684]), nozzleExhaustChokeable(A_fixed = 0.23024082596048073), HP_shaft(J = 1), LP_shaft(J = 1), LPC(N_n_design = 1.0));
  equation
  
    annotation(
      Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-100, 100}, {100, -100}}, textString = "OFF
    DES")}),
      experiment(StartTime = 0, StopTime = 250, Tolerance = 1e-06, Interval = 0.5));
  end OffDesign_Case1_new;
  annotation(
    Icon(graphics = {Text(origin = {3, 4}, extent = {{-63, 54}, {63, -54}}, textString = "J57")}),
    uses(BasicAeroEngines(version = "2.0.0"), Modelica(version = "4.0.0")));
end Assignment_J57;

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "thermalH2TestApp.h"
#include "thermalH2App.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
thermalH2TestApp::validParams()
{
  InputParameters params = thermalH2App::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

thermalH2TestApp::thermalH2TestApp(InputParameters parameters) : MooseApp(parameters)
{
  thermalH2TestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

thermalH2TestApp::~thermalH2TestApp() {}

void
thermalH2TestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  thermalH2App::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"thermalH2TestApp"});
    Registry::registerActionsTo(af, {"thermalH2TestApp"});
  }
}

void
thermalH2TestApp::registerApps()
{
  registerApp(thermalH2App);
  registerApp(thermalH2TestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
thermalH2TestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  thermalH2TestApp::registerAll(f, af, s);
}
extern "C" void
thermalH2TestApp__registerApps()
{
  thermalH2TestApp::registerApps();
}

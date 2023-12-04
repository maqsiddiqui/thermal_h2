#include "thermalH2App.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
thermalH2App::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

thermalH2App::thermalH2App(InputParameters parameters) : MooseApp(parameters)
{
  thermalH2App::registerAll(_factory, _action_factory, _syntax);
}

thermalH2App::~thermalH2App() {}

void 
thermalH2App::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<thermalH2App>(f, af, s);
  Registry::registerObjectsTo(f, {"thermalH2App"});
  Registry::registerActionsTo(af, {"thermalH2App"});

  /* register custom execute flags, action syntax, etc. here */
}

void
thermalH2App::registerApps()
{
  registerApp(thermalH2App);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
thermalH2App__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  thermalH2App::registerAll(f, af, s);
}
extern "C" void
thermalH2App__registerApps()
{
  thermalH2App::registerApps();
}

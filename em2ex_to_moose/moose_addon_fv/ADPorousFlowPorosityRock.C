#include "ADPorousFlowPorosityRock.h"

registerMooseObject("PorousFlowApp", ADPorousFlowPorosityRock);

InputParameters
ADPorousFlowPorosityRock::validParams()
{
  InputParameters params = ADPorousFlowPorosityBase::validParams();
  params.addRequiredCoupledVar("porosity",
                               "Porosity at the reference pressure (a number or a constant "
                               "monomial / finite-volume variable)");
  params.addRequiredParam<Real>("reference_pressure", "Eclipse ROCK reference pressure [Pa]");
  params.addRequiredRangeCheckedParam<Real>(
      "compressibility", "compressibility >= 0", "Eclipse ROCK rock compressibility [1/Pa]");
  params.addClassDescription("Porosity phi_ref (1 + X + X^2/2), X = c_r (P - p_ref), as Eclipse's ROCK");
  return params;
}

ADPorousFlowPorosityRock::ADPorousFlowPorosityRock(const InputParameters & parameters)
  : ADPorousFlowPorosityBase(parameters),
    _phi_ref(coupledValue("porosity")),
    _p_ref(getParam<Real>("reference_pressure")),
    _c_r(getParam<Real>("compressibility")),
    _pf(getADMaterialProperty<Real>(_nodal_material ? "PorousFlow_effective_fluid_pressure_nodal"
                                                    : "PorousFlow_effective_fluid_pressure_qp"))
{
}

void
ADPorousFlowPorosityRock::initQpStatefulProperties()
{
  // the starting ("old") porosity must be the one at the initial pressure, or the first time step
  // sees a spurious pore-volume change
  computeQpProperties();
}

void
ADPorousFlowPorosityRock::computeQpProperties()
{
  const ADReal x = _c_r * (_pf[_qp] - _p_ref);
  _porosity[_qp] = _phi_ref[0] * (1.0 + x + 0.5 * x * x);
}

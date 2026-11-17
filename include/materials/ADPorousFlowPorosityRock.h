#pragma once

#include "PorousFlowPorosityBase.h"

/**
 * Porosity with Eclipse's ROCK compressibility:
 *   phi = phi_ref * (1 + X + X^2 / 2),   X = c_r * (P - p_ref)
 * where P is the effective fluid pressure (saturation-weighted phase pressure, from
 * ADPorousFlowEffectiveFluidPressure), phi_ref the porosity at p_ref (the deck porosity) and
 * c_r the rock compressibility. AD only (cell-centred finite-volume PorousFlow).
 */
class ADPorousFlowPorosityRock : public ADPorousFlowPorosityBase
{
public:
  static InputParameters validParams();
  ADPorousFlowPorosityRock(const InputParameters & parameters);

protected:
  void initQpStatefulProperties() override;
  void computeQpProperties() override;

  const VariableValue & _phi_ref;
  const Real _p_ref;
  const Real _c_r;
  const ADMaterialProperty<Real> & _pf;
};

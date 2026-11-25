#pragma once

#include "FVElementalKernel.h"
#include <unordered_map>

/**
 * Advective flux of one fluid component through non-neighbour connections (NNCs), such as the
 * connections Eclipse makes between cells on the two sides of a fault with throw.
 *
 * Residual of cell a (flux out of a, per unit volume of a):
 *   1/V_a * sum_b T_ab * sum_p lambda_p,up * (Phi_p,a - Phi_p,b)
 *   Phi_p,a - Phi_p,b = (p_a + o_p,a) - (p_b + o_p,b) + rho_p,avg * g . (x_b - x_a)
 *   lambda_p = X_p,c * kr_p * rho_p / mu_p   (upwind)
 * T_ab [m^3] is the Eclipse transmissibility in SI (no viscosity), p the primary pressure
 * (evaluated implicitly at both cells) and o_p a per-phase offset (e.g. capillary pressure).
 * The phase properties (kr, rho, mu, X, o) are functors, typically lagged auxiliary variables.
 */
class FVPorousFlowNNCFlux : public FVElementalKernel
{
public:
  static InputParameters validParams();
  FVPorousFlowNNCFlux(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  const Moose::Functor<ADReal> & _p;
  std::vector<const Moose::Functor<ADReal> *> _offset;
  std::vector<const Moose::Functor<ADReal> *> _kr;
  std::vector<const Moose::Functor<ADReal> *> _rho;
  std::vector<const Moose::Functor<ADReal> *> _mu;
  std::vector<const Moose::Functor<ADReal> *> _x;
  const RealVectorValue _gravity;
  const Real _multiplier;
  std::unordered_map<dof_id_type, std::vector<std::pair<dof_id_type, Real>>> _nnc;
};

#pragma once

#include "FVPorousFlowAdvectiveFlux.h"

/**
 * Two-point (Eclipse-style) advective flux for cell-centred PorousFlow.
 *
 * On internal faces:
 *   flux out of the element = sum over phases of
 *       X * kr * rho / mu |upwind  *  T * Phi
 *   T   = 1 / (1/tE + 1/tN),  t = k_dir (n . D) / (D . D)   (Eclipse's corner-point transmissibility)
 *   Phi = (pE - pN) + rho_face * g . (xN - xE)
 * where D runs from the cell centre to the face centre and k_dir is the permeability along the grid
 * direction of the connection (vertical permeability for mostly-horizontal faces, otherwise the
 * horizontal permeability along D). On Croft this reproduces Eclipse's TRANX/Y/Z. Gravity is taken along the
 * same centre-to-centre line as the pressure difference, so a hydrostatic state gives zero
 * flux on any grid geometry, including dipping and distorted corner-point cells.
 * Boundary faces are treated as in FVPorousFlowAdvectiveFlux.
 */
class FVPorousFlowAdvectiveFluxTPFA : public FVPorousFlowAdvectiveFlux
{
public:
  static InputParameters validParams();
  FVPorousFlowAdvectiveFluxTPFA(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;
};

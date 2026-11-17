#include "FVPorousFlowAdvectiveFluxTPFA.h"
#include "PorousFlowDictator.h"

registerADMooseObject("PorousFlowApp", FVPorousFlowAdvectiveFluxTPFA);

InputParameters
FVPorousFlowAdvectiveFluxTPFA::validParams()
{
  InputParameters params = FVPorousFlowAdvectiveFlux::validParams();
  params.addClassDescription("Two-point (harmonic transmissibility, hydrostatically consistent) "
                             "advective Darcy flux for cell-centred PorousFlow");
  return params;
}

FVPorousFlowAdvectiveFluxTPFA::FVPorousFlowAdvectiveFluxTPFA(const InputParameters & params)
  : FVPorousFlowAdvectiveFlux(params)
{
}

ADReal
FVPorousFlowAdvectiveFluxTPFA::computeQpResidual()
{
  if (onBoundary(*_face_info))
    return FVPorousFlowAdvectiveFlux::computeQpResidual();

  const auto & fi = *_face_info;
  const RealVectorValue n = fi.normal(); // unit normal, element -> neighbor
  const RealVectorValue xE = fi.elemCentroid();
  const RealVectorValue xN = fi.neighborCentroid();
  const RealVectorValue xf = fi.faceCentroid();
  const RealVectorValue DE = xf - xE;
  const RealVectorValue DN = xN - xf;

  // Eclipse's corner-point transmissibility: T = 1 / (1/tE + 1/tN), t = k_dir (n . D) / (D . D) per unit
  // area, with k_dir the permeability along the grid direction of the connection: the vertical
  // permeability for faces that are mostly horizontal (K connections), otherwise the horizontal
  // permeability along the horizontal direction from the cell centre to the face (I/J connections).
  const bool vertical = std::abs(n(2)) > 0.70710678;
  auto kdir = [&](const auto & K, const RealVectorValue & D)
  {
    if (vertical)
      return K(2, 2);
    RealVectorValue h(D(0), D(1), 0.0);
    const Real hn = h.norm();
    if (hn <= 0.0)
      return K(2, 2);
    h /= hn;
    return h * (K * h);
  };
  const ADReal kE = kdir(_permeability[_qp], DE);
  const ADReal kN = kdir(_permeability_neighbor[_qp], DN);
  const Real nDE = n * DE, nDN = n * DN;
  if (kE <= 0.0 || kN <= 0.0 || nDE <= 0.0 || nDN <= 0.0)
    return 0.0;
  const ADReal tE = kE * nDE / (DE * DE);
  const ADReal tN = kN * nDN / (DN * DN);
  const ADReal trans = 1.0 / (1.0 / tE + 1.0 / tN); // m^2 / m

  const Real g_dx = _gravity * (xN - xE);

  ADReal flux = 0.0;
  for (const auto p : make_range(_num_phases))
  {
    const ADReal rho_face = 0.5 * (_density[_qp][p] + _density_neighbor[_qp][p]);
    const ADReal phi = (_pressure[_qp][p] - _pressure_neighbor[_qp][p]) + rho_face * g_dx;
    ADReal lambda;
    if (phi > 0.0)
      lambda = _mass_fractions[_qp][p][_fluid_component] * _relperm[_qp][p] * _density[_qp][p] /
               _viscosity[_qp][p];
    else
      lambda = _mass_fractions_neighbor[_qp][p][_fluid_component] * _relperm_neighbor[_qp][p] *
               _density_neighbor[_qp][p] / _viscosity_neighbor[_qp][p];
    flux += lambda * trans * phi;
  }
  return flux;
}

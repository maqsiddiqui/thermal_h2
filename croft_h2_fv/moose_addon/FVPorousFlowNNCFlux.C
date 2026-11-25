#include "FVPorousFlowNNCFlux.h"
#include "NNCRelationshipManager.h"
#include "Assembly.h"
#include "MooseMesh.h"

registerADMooseObject("PorousFlowApp", FVPorousFlowNNCFlux);

InputParameters
FVPorousFlowNNCFlux::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription("Advective flux of a fluid component through non-neighbour "
                             "connections (two-point, upwind, transmissibilities from a file)");
  params.addRequiredParam<FileName>("nnc_file", "CSV file: elem_a, elem_b, trans [m^3]");
  params.addRequiredParam<MooseFunctorName>("pressure", "Primary pressure (implicit)");
  params.addRequiredParam<std::vector<MooseFunctorName>>(
      "pressure_offsets", "Per phase: phase pressure minus the primary pressure (e.g. 0, Pc)");
  params.addRequiredParam<std::vector<MooseFunctorName>>("relperm", "Per phase relative permeability");
  params.addRequiredParam<std::vector<MooseFunctorName>>("density", "Per phase density");
  params.addRequiredParam<std::vector<MooseFunctorName>>("viscosity", "Per phase viscosity");
  params.addRequiredParam<std::vector<MooseFunctorName>>(
      "mass_fraction", "Per phase mass fraction of this component");
  params.addRequiredParam<RealVectorValue>("gravity", "Gravitational acceleration");
  params.addParam<Real>("multiplier", 1.0, "Multiplier on all NNC transmissibilities");
  params.addRelationshipManager(
      "NNCRelationshipManager",
      Moose::RelationshipManagerType::ALGEBRAIC | Moose::RelationshipManagerType::COUPLING,
      [](const InputParameters & obj_params, InputParameters & rm_params)
      { rm_params.set<FileName>("nnc_file") = obj_params.get<FileName>("nnc_file"); });
  return params;
}

FVPorousFlowNNCFlux::FVPorousFlowNNCFlux(const InputParameters & params)
  : FVElementalKernel(params),
    _p(getFunctor<ADReal>("pressure")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _multiplier(getParam<Real>("multiplier"))
{
  auto grab = [&](const std::string & name, std::vector<const Moose::Functor<ADReal> *> & v)
  {
    for (const auto & f : getParam<std::vector<MooseFunctorName>>(name))
      v.push_back(&getFunctor<ADReal>(f));
  };
  grab("pressure_offsets", _offset);
  grab("relperm", _kr);
  grab("density", _rho);
  grab("viscosity", _mu);
  grab("mass_fraction", _x);
  const auto np = _offset.size();
  if (_kr.size() != np || _rho.size() != np || _mu.size() != np || _x.size() != np)
    paramError("relperm", "pressure_offsets, relperm, density, viscosity and mass_fraction must "
                          "all have one entry per phase");
  NNCRelationshipManager::readNNCFile(getParam<FileName>("nnc_file"), _nnc);
}

ADReal
FVPorousFlowNNCFlux::computeQpResidual()
{
  const auto it = _nnc.find(_current_elem->id());
  if (it == _nnc.end())
    return 0.0;

  const auto state = determineState();
  const auto ea = makeElemArg(_current_elem);
  const RealVectorValue xa = _current_elem->vertex_average();
  const ADReal pa = _p(ea, state);

  ADReal flux = 0.0;
  for (const auto & [bid, trans] : it->second)
  {
    const Elem * const elem_b = _mesh.getMesh().query_elem_ptr(bid);
    if (!elem_b)
      mooseError("FVPorousFlowNNCFlux: element ", bid, " is not available on this processor");
    const auto eb = makeElemArg(elem_b);
    const Real g_dx = _gravity * (elem_b->vertex_average() - xa);
    const ADReal pb = _p(eb, state);
    for (const auto p : index_range(_kr))
    {
      const ADReal rho_a = (*_rho[p])(ea, state), rho_b = (*_rho[p])(eb, state);
      const ADReal phi = (pa + (*_offset[p])(ea, state)) - (pb + (*_offset[p])(eb, state)) +
                         0.5 * (rho_a + rho_b) * g_dx;
      ADReal lambda;
      if (phi > 0.0)
        lambda = (*_x[p])(ea, state) * (*_kr[p])(ea, state) * rho_a / (*_mu[p])(ea, state);
      else
        lambda = (*_x[p])(eb, state) * (*_kr[p])(eb, state) * rho_b / (*_mu[p])(eb, state);
      flux += lambda * _multiplier * trans * phi;
    }
  }
  return flux / _assembly.elemVolume();
}

//* Add-on for the Croft H2 storage model (finite volume). LGPL 2.1, as MOOSE.
#include "FVMaxIncrementDamper.h"
#include "SystemBase.h"
#include "MooseMesh.h"
#include "MooseVariableFieldBase.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"

registerMooseObject("PorousFlowApp", FVMaxIncrementDamper);

InputParameters
FVMaxIncrementDamper::validParams()
{
  InputParameters params = GeneralDamper::validParams();
  params.addClassDescription("Limits the Newton change of an elemental (e.g. finite-volume) variable "
                             "per iteration, and optionally keeps it within bounds.");
  params.addRequiredParam<NonlinearVariableName>("variable", "The variable to damp");
  params.addParam<Real>("max_increment", std::numeric_limits<Real>::max(),
                        "Largest allowed change of the variable in one Newton iteration");
  params.addParam<Real>("min_value", std::numeric_limits<Real>::lowest(), "Lower bound");
  params.addParam<Real>("max_value", std::numeric_limits<Real>::max(), "Upper bound");
  params.addParam<NonlinearVariableName>(
      "ignore_variable", "Skip cells where this variable is below ignore_below (e.g. a composition "
      "variable in cells with almost no gas, where it is poorly determined and its Newton changes are "
      "large but harmless)");
  params.addParam<Real>("ignore_below", 0.0, "Threshold for ignore_variable");
  params.addParam<bool>("verbose", false, "Print the cell that limits the damping");
  params.addRangeCheckedParam<Real>("min_bound_damping", 0.1, "min_bound_damping > 0 & min_bound_damping <= 1",
                                    "The bounds never damp the update below this factor");
  return params;
}

FVMaxIncrementDamper::FVMaxIncrementDamper(const InputParameters & parameters)
  : GeneralDamper(parameters),
    _var_num(_sys.getVariable(0, getParam<NonlinearVariableName>("variable")).number()),
    _max_increment(getParam<Real>("max_increment")),
    _min_value(getParam<Real>("min_value")),
    _max_value(getParam<Real>("max_value")),
    _min_bound_damping(getParam<Real>("min_bound_damping")),
    _ignore_var_num(isParamValid("ignore_variable")
                        ? int(_sys.getVariable(0, getParam<NonlinearVariableName>("ignore_variable")).number())
                        : -1),
    _ignore_below(getParam<Real>("ignore_below"))
{
}

Real
FVMaxIncrementDamper::computeDamping(const NumericVector<Number> & solution,
                                     const NumericVector<Number> & update)
{
  // MOOSE: new = old - damping * update, and 'solution' is the undamped new iterate
  Real damping = 1.0;
  Real worst_u = 0, worst_du = 0;
  Point worst_x;
  const libMesh::DofMap & dof_map = _sys.dofMap();
  std::vector<libMesh::dof_id_type> dofs, idofs;
  const auto first = solution.first_local_index(), last = solution.last_local_index();
  for (const auto & elem : *_sys.mesh().getActiveLocalElementRange())
  {
    dof_map.dof_indices(elem, dofs, _var_num);
    if (_ignore_var_num >= 0)
    {
      dof_map.dof_indices(elem, idofs, (unsigned int)_ignore_var_num);
      bool skip = false;
      for (const auto j : idofs)
        if (j >= first && j < last && solution(j) + update(j) < _ignore_below) // old value of that variable
          skip = true;
      if (skip)
        continue;
    }
    for (const auto i : dofs)
    {
      if (i < first || i >= last)
        continue;
      const Real du = -update(i);              // full Newton change
      const Real u_old = solution(i) - du;
      const Real adu = std::abs(du);
      const Real d_before = damping;
      if (adu > _max_increment)
        damping = std::min(damping, _max_increment / adu);
      const Real u_new = u_old + damping * du;
      if (u_old >= _min_value && u_new < _min_value && du < 0)
        damping = std::min(damping, std::max(_min_bound_damping, (_min_value - u_old) / du));
      if (u_old <= _max_value && u_new > _max_value && du > 0)
        damping = std::min(damping, std::max(_min_bound_damping, (_max_value - u_old) / du));
      if (damping < d_before)
      {
        worst_u = u_old;
        worst_du = du;
        worst_x = elem->vertex_average();
      }
    }
  }
  if (getParam<bool>("verbose") && damping < 0.5)
    std::cout << "  [rank " << processor_id() << "] " << name() << ": damping " << damping << " set by the cell at " << worst_x
               << " (value " << worst_u << ", Newton change " << worst_du << ")" << std::endl;
  _communicator.min(damping);
  return damping;
}

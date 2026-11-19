//* Add-on for the Croft H2 storage model (finite volume). LGPL 2.1, as MOOSE.
#pragma once

#include "GeneralDamper.h"

/**
 * Newton damping for finite-volume (or any elemental) variables, which MOOSE's MaxIncrement
 * (an ElementDamper for finite-element variables) does not handle. The whole Newton update is scaled
 * by one factor (the smallest over all cells and processors) so that
 *   |change of the variable| <= max_increment in every cell, and
 *   a cell inside [min_value, max_value] is not pushed further than to that bound
 *   (the bound part never damps below min_bound_damping, so Newton cannot stall on it).
 */
class FVMaxIncrementDamper : public GeneralDamper
{
public:
  static InputParameters validParams();
  FVMaxIncrementDamper(const InputParameters & parameters);

  virtual Real computeDamping(const NumericVector<Number> & solution,
                              const NumericVector<Number> & update) override;

protected:
  const unsigned int _var_num;
  const Real _max_increment;
  const Real _min_value;
  const Real _max_value;
  const Real _min_bound_damping;
  const int _ignore_var_num;
  const Real _ignore_below;
};

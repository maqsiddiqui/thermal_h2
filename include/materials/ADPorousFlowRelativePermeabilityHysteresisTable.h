//* Add-on for the Croft models (finite volume). LGPL 2.1, as MOOSE.
#pragma once

#include "PorousFlowRelativePermeabilityBase.h"

/**
 * Tabulated relative permeability of one phase with Eclipse's two-phase hysteresis options
 * (SATOPTS HYSTER, EHYSTR item 2), as in ECLIPSE and OPM Flow:
 *
 *   hysteresis = none        kr = m * kr_drainage(S)
 *   hysteresis = imbibition  kr = m * kr_imbibition(S)      (wetting phase, EHYSTR options 1 and 3)
 *   hysteresis = killough    non-wetting phase, Killough (1976) scanning curves (EHYSTR 2 and 3):
 *        S >= S_hy:  kr = m * kr_d(S)                         (on the drainage curve)
 *        S <  S_hy:  S_t   = S_dc + (S_hy - S_dc) / (1 + C (S_hy - S_dc))          (Land trapping)
 *                    C     = 1 / (S_ic - S_dc) - 1 / (S_max - S_dc)
 *                    S_n   = S_ic + (S - S_t) (S_max - S_ic) / (S_hy - S_t)
 *                    kr    = m * kr_i(S_n) * kr_d(S_hy) / kr_i(S_max)
 *        S_dc, S_ic: critical saturations of the drainage and imbibition tables (largest S with kr = 0),
 *        S_max: largest saturation of the drainage table, S_hy: largest saturation the cell has had
 *        (coupled variable max_saturation, updated once per time step, i.e. lagged).
 *
 * m (kr_multiplier) is a per-cell multiplier, e.g. Eclipse's IKRW / KRW end-point scaling.
 * Tables are interpolated linearly and held constant beyond their ends, as in Eclipse, except that every
 * corner is rounded: interior points and the last point by a quadratic over +-smoothing (default 0.005,
 * at most a quarter of the neighbouring intervals), the first point (critical saturation, below which kr
 * stays 0) by a cubic over [x0, x0 + 2 smoothing], and the switch from the scanning curve to the drainage
 * curve at S_hy is blended over switch_width: both make kr continuously differentiable, which Newton needs
 * (with sharp corners it can cycle between the two sides of a corner). smoothing = switch_width = 0 gives
 * Eclipse's curves exactly.
 * AD only.
 */
class ADPorousFlowRelativePermeabilityHysteresisTable : public ADPorousFlowRelativePermeabilityBase
{
public:
  static InputParameters validParams();
  ADPorousFlowRelativePermeabilityHysteresisTable(const InputParameters & parameters);

protected:
  void computeQpProperties() override;
  ADReal relativePermeability(ADReal) const override { return 0.0; } // not used
  Real dRelativePermeability(Real) const override { return 0.0; }    // not used

  /// linear interpolation, constant outside the table, corners rounded over +-h
  template <typename T>
  static T interp(const std::vector<Real> & x, const std::vector<Real> & y, const T & s, Real h);
  /// largest saturation with kr = 0 (the critical saturation), or the first saturation
  static Real critical(const std::vector<Real> & x, const std::vector<Real> & y);
  static void check(const std::vector<Real> & x, const std::vector<Real> & y, const std::string & what,
                    const MooseObject & obj);

  enum class Model { NONE, IMBIBITION, KILLOUGH } _model;
  const std::vector<Real> _sd, _krd; ///< drainage table
  const std::vector<Real> _si, _kri; ///< imbibition table
  const VariableValue & _s_hy;       ///< historical maximum saturation of this phase
  const VariableValue & _mult;       ///< per-cell multiplier
  const Real _h;                     ///< corner-rounding half width
  const Real _switch;                ///< scanning -> drainage blending width
  Real _s_dc = 0, _s_ic = 0, _s_max = 1, _C = 0, _kri_max = 1;
};

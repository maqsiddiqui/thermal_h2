//* Add-on for the Croft models (finite volume). LGPL 2.1, as MOOSE.
#include "ADPorousFlowRelativePermeabilityHysteresisTable.h"

registerMooseObject("PorousFlowApp", ADPorousFlowRelativePermeabilityHysteresisTable);

InputParameters
ADPorousFlowRelativePermeabilityHysteresisTable::validParams()
{
  InputParameters params = ADPorousFlowRelativePermeabilityBase::validParams();
  params.addRequiredParam<std::vector<Real>>("saturation",
                                             "Drainage table: saturation of this phase (increasing)");
  params.addRequiredParam<std::vector<Real>>("kr", "Drainage table: relative permeability");
  params.addParam<std::vector<Real>>("imbibition_saturation", {},
                                     "Imbibition table: saturation of this phase (increasing)");
  params.addParam<std::vector<Real>>("imbibition_kr", {}, "Imbibition table: relative permeability");
  params.addParam<MooseEnum>("hysteresis", MooseEnum("none imbibition killough", "none"),
                             "none: drainage table; imbibition: imbibition table (Eclipse EHYSTR 1, 3 "
                             "for the wetting phase); killough: Killough scanning curves (EHYSTR 2, 3 "
                             "for the non-wetting phase)");
  params.addCoupledVar("max_saturation", 0.0,
                       "Largest saturation of this phase the cell has had (killough); update it once per "
                       "time step, e.g. a ParsedAux max(itself, saturation) executed on timestep_end");
  params.addCoupledVar("kr_multiplier", 1.0, "Per-cell multiplier of the relative permeability");
  params.addRangeCheckedParam<Real>("smoothing", 0.005, "smoothing >= 0",
                                    "Half width (in saturation) over which the corners of the tables are rounded "
                                    "(limited to a quarter of the neighbouring intervals); 0: linear, as Eclipse");
  params.addRangeCheckedParam<Real>("switch_width", 0.01, "switch_width >= 0",
                                    "killough: width (in saturation) below S_hy over which the scanning curve "
                                    "is blended into the drainage curve; 0: sharp switch, as Eclipse");
  params.addClassDescription("Tabulated relative permeability with Eclipse's hysteresis options "
                             "(imbibition curve, Killough scanning curves)");
  return params;
}

ADPorousFlowRelativePermeabilityHysteresisTable::ADPorousFlowRelativePermeabilityHysteresisTable(
    const InputParameters & parameters)
  : ADPorousFlowRelativePermeabilityBase(parameters),
    _model(Model(int(getParam<MooseEnum>("hysteresis")))),
    _sd(getParam<std::vector<Real>>("saturation")),
    _krd(getParam<std::vector<Real>>("kr")),
    _si(getParam<std::vector<Real>>("imbibition_saturation")),
    _kri(getParam<std::vector<Real>>("imbibition_kr")),
    _s_hy(coupledValue("max_saturation")),
    _mult(coupledValue("kr_multiplier")),
    _h(getParam<Real>("smoothing")),
    _switch(getParam<Real>("switch_width"))
{
  check(_sd, _krd, "saturation / kr", *this);
  if (_model != Model::NONE)
    check(_si, _kri, "imbibition_saturation / imbibition_kr", *this);
  if (_model == Model::KILLOUGH && !isCoupled("max_saturation"))
    paramError("max_saturation", "is needed for hysteresis = killough");
  if (_model == Model::KILLOUGH)
  {
    _s_dc = critical(_sd, _krd);
    _s_ic = critical(_si, _kri);
    _s_max = _sd.back();
    _kri_max = interp<Real>(_si, _kri, _s_max, _h);
    if (_s_ic > _s_dc && _s_max > _s_ic && _kri_max > 0)
      _C = 1.0 / (_s_ic - _s_dc) - 1.0 / (_s_max - _s_dc);
    else
      mooseWarning(name(), ": the imbibition table traps no more than the drainage table; "
                           "Killough hysteresis has no effect");
  }
}

void
ADPorousFlowRelativePermeabilityHysteresisTable::check(const std::vector<Real> & x,
                                                       const std::vector<Real> & y,
                                                       const std::string & what,
                                                       const MooseObject & obj)
{
  if (x.size() < 2 || x.size() != y.size())
    obj.mooseError(what, ": need two equally long lists with at least 2 entries");
  for (std::size_t i = 1; i < x.size(); ++i)
    if (!(x[i] > x[i - 1]))
      obj.mooseError(what, ": saturations must increase");
}

Real
ADPorousFlowRelativePermeabilityHysteresisTable::critical(const std::vector<Real> & x,
                                                          const std::vector<Real> & y)
{
  Real sc = x.front();
  for (std::size_t i = 0; i < x.size(); ++i)
    if (y[i] <= 0.0)
      sc = x[i];
    else
      break;
  return sc;
}

template <typename T>
T
ADPorousFlowRelativePermeabilityHysteresisTable::interp(const std::vector<Real> & x,
                                                        const std::vector<Real> & y,
                                                        const T & s,
                                                        Real h)
{
  const std::size_t n = x.size();
  const Real sv = MetaPhysicL::raw_value(s);
  // slope left and right of table point k (0 outside the table)
  auto slope = [&](std::size_t i) { return (y[i] - y[i - 1]) / (x[i] - x[i - 1]); }; // segment i-1 .. i
  if (h > 0.0)
  {
    // the corner whose rounding window contains s, if any
    std::size_t k = std::lower_bound(x.begin(), x.end(), sv) - x.begin(); // x[k-1] < s <= x[k]
    for (std::size_t c : {k == 0 ? std::size_t(0) : k - 1, std::min(k, n - 1)})
    {
      if (c == 0)
      {
        // first point (critical saturation): kr stays y[0] below it; above it a cubic leaves with zero
        // slope and joins the first segment with its slope at x[0] + 2 hc
        const Real hc = std::min(h, 0.25 * (x[1] - x[0]));
        const T d = s - x[0];
        if (sv >= x[0] && sv < x[0] + 2.0 * hc)
        {
          const Real m = slope(1);
          return y[0] + m * d * d / hc - m * d * d * d / (4.0 * hc * hc);
        }
        continue;
      }
      const Real wl = c == 0 ? x[1] - x[0] : x[c] - x[c - 1];
      const Real wr = c == n - 1 ? x[n - 1] - x[n - 2] : x[c + 1] - x[c];
      const Real hc = std::min(h, 0.25 * std::min(wl, wr));
      if (std::abs(sv - x[c]) < hc)
      {
        const Real ml = c == 0 ? 0.0 : slope(c);
        const Real mr = c == n - 1 ? 0.0 : slope(c + 1);
        const T d = s - x[c];
        return y[c] + ml * d + (mr - ml) * (d + hc) * (d + hc) / (4.0 * hc);
      }
    }
  }
  if (sv <= x.front())
    return y.front();
  if (sv >= x.back())
    return y.back();
  const auto i = std::upper_bound(x.begin(), x.end(), sv) - x.begin(); // x[i-1] <= s < x[i]
  return y[i - 1] + (y[i] - y[i - 1]) * (s - x[i - 1]) / (x[i] - x[i - 1]);
}

void
ADPorousFlowRelativePermeabilityHysteresisTable::computeQpProperties()
{
  const ADReal & s = _saturation[_qp][_phase_num];
  ADReal kr;
  switch (_model)
  {
    case Model::NONE:
      kr = interp(_sd, _krd, s, _h);
      break;
    case Model::IMBIBITION:
      kr = interp(_si, _kri, s, _h);
      break;
    case Model::KILLOUGH:
    {
      const Real s_hy = std::min(_s_hy[_qp], _s_max);
      const Real sv = MetaPhysicL::raw_value(s);
      const ADReal kr_d = interp(_sd, _krd, s, _h);
      if (_C <= 0.0 || s_hy <= _s_dc || sv >= s_hy)
        kr = kr_d; // on the drainage curve
      else
      {
        const Real s_t = _s_dc + (s_hy - _s_dc) / (1.0 + _C * (s_hy - _s_dc)); // trapped saturation
        if (s_hy - s_t <= 1e-9)
          kr = 0.0;
        else
        {
          const ADReal s_n = _s_ic + (s - s_t) * (_s_max - _s_ic) / (s_hy - s_t);
          const ADReal kr_s = interp(_si, _kri, s_n, _h) * (interp<Real>(_sd, _krd, s_hy, _h) / _kri_max);
          const Real w0 = std::min(_switch, 0.5 * (s_hy - s_t)); // blend into the drainage curve below S_hy
          if (w0 > 0.0 && sv > s_hy - w0)
          {
            const ADReal u = (s - (s_hy - w0)) / w0;
            const ADReal w = u * u * (3.0 - 2.0 * u);
            kr = (1.0 - w) * kr_s + w * kr_d;
          }
          else
            kr = kr_s;
        }
      }
      break;
    }
  }
  _relative_permeability[_qp] = kr * _mult[_qp] * _scaling;
}

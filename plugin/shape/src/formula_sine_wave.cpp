#include <cmath>
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/formula_sine_wave.h"

namespace feasst {

FormulaSineWave::FormulaSineWave(argtype args) : Formula(&args) {
  class_name_ = "FormulaSineWave";
  amplitude_ = dble("amplitude", &args, 1.);
  width_ = dble("width", &args, 2*PI);
  phase_ = dble("phase", &args, 0.);
  shift_ = dble("shift", &args, default_shift());
  FEASST_CHECK_ALL_USED(args);
}

class MapFormulaSineWave {
 public:
  MapFormulaSineWave() {
    FormulaSineWave().deserialize_map()["FormulaSineWave"] =
      std::make_shared<FormulaSineWave>();
  }
};

static MapFormulaSineWave mapper_ = MapFormulaSineWave();

std::shared_ptr<Formula> FormulaSineWave::create(std::istream& istr) const {
  return std::make_shared<FormulaSineWave>(istr);
}

FormulaSineWave::FormulaSineWave(std::istream& istr)
  : Formula(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5404, "version mismatch: " << version);
  feasst_deserialize(&amplitude_, istr);
  feasst_deserialize(&width_, istr);
  feasst_deserialize(&phase_, istr);
  feasst_deserialize(&shift_, istr);
}

void FormulaSineWave::serialize(std::ostream& ostr) const {
  serialize_formula_(ostr);
  feasst_serialize_version(5404, ostr);
  feasst_serialize(amplitude_, ostr);
  feasst_serialize(width_, ostr);
  feasst_serialize(phase_, ostr);
  feasst_serialize(shift_, ostr);
}

double FormulaSineWave::evaluate(const double x) const {
  return shift_ + amplitude_*std::sin(2*PI*(x - phase_)/width_);
}

double FormulaSineWave::derivative(const double x) const {
  return amplitude_*2*PI/width_*std::cos(2*PI*(x - phase_)/width_);
}

double FormulaSineWave::second_derivative(const double x) const {
  return -amplitude_*std::pow(2*PI/width_, 2)*std::sin(2*PI*(x - phase_)/width_);
}

double FormulaSineWave::nearest_minimum(const double x) const {
  TRACE("x " << x);
  TRACE(x-phase_-0.75*width_);
  TRACE((x-phase_-0.75*width_)/width_);
  const int iter = feasst::round((x - phase_ - 0.75*width_)/width_);
  TRACE("iter " << iter);
  return 0.75*width_ + iter*width_ + phase_;
}

double FormulaSineWave::nearest_maximum(const double x) const {
  TRACE("x " << x);
  TRACE(x-phase_-0.75*width_);
  TRACE((x-phase_-0.75*width_)/width_);
  const int iter = feasst::round((x - phase_ - 0.25*width_)/width_);
  TRACE("iter " << iter);
  return 0.25*width_ + iter*width_ + phase_;
}

}  // namespace feasst

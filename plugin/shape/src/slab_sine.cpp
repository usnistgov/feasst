#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/half_space_sine.h"
#include "shape/include/slab_sine.h"

namespace feasst {

FEASST_MAPPER(SlabSine, argtype({{"dimension", "0"}, {"wave_dimension", "0"},
  {"average_bound0", "-1"}, {"average_bound1", "1"}}));

SlabSine::SlabSine(argtype * args) : ShapeIntersect() {
  class_name_ = "SlabSine";
  int dimension = integer("dimension", args);
  int wave_dimension = integer("wave_dimension", args);
  double upper = dble("average_bound0", args);
  double lower = dble("average_bound1", args);
  feasst_sort(&lower, &upper);
  ASSERT(upper - lower > NEAR_ZERO, "slab is too thin");
  auto sine_wave = std::make_shared<FormulaSineWave>(args);
  sine_wave->set_phase(sine_wave->phase() - 0.25*sine_wave->width());
  auto half0 = MakeHalfSpaceSine(sine_wave, {
    {"dimension", str(dimension)},
    {"wave_dimension", str(wave_dimension)},
    {"intersection", str(upper)},
    {"direction", str(-1)},
  });
  sine_wave->set_phase(sine_wave->phase() + 0.5*sine_wave->width());
  auto half1 = MakeHalfSpaceSine(sine_wave, {
    {"dimension", str(dimension)},
    {"wave_dimension", str(wave_dimension)},
    {"intersection", str(lower)},
    {"direction", str(1)},
  });
  set(half0, half1);
}
SlabSine::SlabSine(argtype args) : SlabSine(&args) {
  feasst_check_all_used(args);
}

void SlabSine::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_shape_intersect_(ostr);
  feasst_serialize_version(7689, ostr);
}

SlabSine::SlabSine(std::istream& istr) : ShapeIntersect(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(7689 == version, version);
}

}  // namespace feasst

#include "utils/include/utils_io.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/half_space_sine.h"
#include "shape/include/slab_sine.h"

namespace feasst {

SlabSine::SlabSine(std::shared_ptr<FormulaSineWave> sine_wave,
    const argtype &args) : ShapeIntersect() {
  Arguments args_(args);
  int dimension = args_.key("dimension").integer();
  int wave_dimension = args_.key("wave_dimension").integer();
  double upper = args_.key("average_bound0").dble();
  double lower = args_.key("average_bound1").dble();
  sort(&lower, &upper);
  ASSERT(upper - lower > NEAR_ZERO, "slab is too thin");
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

}  // namespace feasst

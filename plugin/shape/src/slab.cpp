#include "utils/include/utils_io.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/half_space.h"
#include "shape/include/slab.h"

namespace feasst {

Slab::Slab(const argtype &args) : ShapeIntersect() {
  Arguments args_(args);
  int dimension = args_.key("dimension").integer();
  double upper = args_.key("bound0").dble();
  double lower = args_.key("bound1").dble();
  sort(&lower, &upper);
  ASSERT(upper - lower > NEAR_ZERO, "slab is infinitesimally thin");
  auto half0 = MakeHalfSpace({
    {"dimension", str(dimension)},
    {"intersection", str(upper)},
    {"direction", str(-1)},
  });
  auto half1 = MakeHalfSpace({
    {"dimension", str(dimension)},
    {"intersection", str(lower)},
    {"direction", str(1)},
  });
  set(half0, half1);
}

}  // namespace feasst

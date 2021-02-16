#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/half_space.h"
#include "shape/include/slab.h"

namespace feasst {

Slab::Slab(argtype args) : ShapeIntersect() {
  int dimension = integer("dimension", &args);
  double upper = dble("bound0", &args);
  double lower = dble("bound1", &args);
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
  check_all_used(args);
}

}  // namespace feasst

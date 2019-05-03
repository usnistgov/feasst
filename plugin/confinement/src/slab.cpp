#include <cmath>
#include "confinement/include/slab.h"
#include "confinement/include/half_space.h"
#include "core/include/utils_math.h"

namespace feasst {

class MapSlab {
 public:
  MapSlab() {
    auto obj = MakeSlab({
      {"dimension", "0"},
      {"bound0", "1"},
      {"bound1", "-1"},
    });
    obj->deserialize_map()["Slab"] = obj;
  }
};

static MapSlab mapper_ = MapSlab();

Slab::Slab(const argtype &args) : Shape() {
  args_.init(args);
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
  slab_ = std::make_shared<ShapeIntersect>(half0, half1);
}

}  // namespace feasst

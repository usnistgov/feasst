#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/half_space.h"
#include "shape/include/slab.h"

namespace feasst {

class MapSlab {
 public:
  MapSlab() {
    auto obj = MakeSlab({{"dimension", "0"}, {"bound0", "0"}, {"bound1", "1"}});
    obj->deserialize_map()["Slab"] = obj;
  }
};

static MapSlab mapper_ = MapSlab();

Slab::Slab(argtype * args) : ShapeIntersect() {
  class_name_ = "Slab";
  int dimension = integer("dimension", args);
  double upper = dble("bound0", args);
  double lower = dble("bound1", args);
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
Slab::Slab(argtype args) : Slab(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void Slab::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_shape_intersect_(ostr);
  feasst_serialize_version(7126, ostr);
}

Slab::Slab(std::istream& istr) : ShapeIntersect(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(7126 == version, version);
}

}  // namespace feasst

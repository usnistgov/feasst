#include <cmath>
#include "utils/include/serialize.h"
#include "confinement/include/shape_intersect.h"
#include "confinement/include/slab.h"
#include "confinement/include/half_space.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"

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

void Slab::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(485, ostr);
  feasst_serialize(slab_, ostr);
}

Slab::Slab(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(485 == version, version);
  // feasst_deserialize(slab_, istr);
  // HWH for unknown reasons, the above doesn't work
  int existing;
  istr >> existing;
  if (existing != 0) {
    slab_ = slab_->deserialize(istr);
  }
}

}  // namespace feasst

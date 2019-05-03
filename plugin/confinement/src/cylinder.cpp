#include <cmath>
#include "confinement/include/cylinder.h"

namespace feasst {

class MapCylinder {
 public:
  MapCylinder() {
    auto obj = MakeCylinder(
      {{"radius", "1"}},
      Position().set_vector({0, 0, 0}),
      Position().set_vector({0, 0, 1})
    );
    obj->deserialize_map()["Cylinder"] = obj;
  }
};

static MapCylinder mapper_ = MapCylinder();

Cylinder::Cylinder(const argtype &args,
    const Position point0,
    const Position point1) : Shape() {
  args_.init(args);
  radius_ = args_.key("radius").dble();
  point0_ = point0;
  point1_ = point1;
}

double Cylinder::nearest_distance(const Position& point) const {
  return point.nearest_distance_to_axis(point0_, point1_) - radius_;
}

}  // namespace feasst

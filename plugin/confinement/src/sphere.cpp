#include <cmath>
#include "confinement/include/sphere.h"

namespace feasst {

class MapSphere {
 public:
  MapSphere() {
    auto obj = MakeSphere(
      {{"radius", "1"}},
      Position().set_vector({0, 0, 0})
    );
    obj->deserialize_map()["Sphere"] = obj;
  }
};

static MapSphere mapper_ = MapSphere();

Sphere::Sphere(const argtype &args,
    const Position center) : Shape() {
  args_.init(args);
  radius_ = args_.key("radius").dble();
  center_ = center;
}

double Sphere::nearest_distance(const Position& point) const {
  Position relative = point;
  relative.subtract(center_);
  return relative.distance() - radius_;
}

}  // namespace feasst

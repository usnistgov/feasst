#include <cmath>
#include "utils/include/serialize.h"
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

void Sphere::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(629, ostr);
  feasst_serialize(radius_, ostr);
  feasst_serialize_fstobj(center_, ostr);
}

Sphere::Sphere(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(629 == version, version);
  feasst_deserialize(&radius_, istr);
  feasst_deserialize_fstobj(&center_, istr);
}

}  // namespace feasst

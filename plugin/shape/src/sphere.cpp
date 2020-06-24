#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/sphere.h"

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
  class_name_ = "Sphere";
  args_.init(args);
  radius_ = args_.key("radius").dble();
  center_ = center;
}

double Sphere::nearest_distance(const Position& point) const {
  return point.distance(center_) - radius_;
}

void Sphere::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_shape_(ostr);
  feasst_serialize_version(629, ostr);
  feasst_serialize(radius_, ostr);
  feasst_serialize_fstobj(center_, ostr);
}

Sphere::Sphere(std::istream& istr) : Shape(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(629 == version, version);
  feasst_deserialize(&radius_, istr);
  feasst_deserialize_fstobj(&center_, istr);
}

double Sphere::surface_area() const { return 4.*PI*std::pow(radius_, 2); }

double Sphere::volume() const { return 4./3.*PI*std::pow(radius_, 3); }

void Sphere::surface_mesh(const int num,
    std::vector<Position> * points) const {
  if (num != static_cast<int>(points->size())) {
    points->resize(num);
    Position origin;
    origin.set_to_origin(3);
    for (int pt = 0; pt < num; ++pt) (*points)[pt] = origin;
  }
  const double increment = PI*(3. - std::sqrt(5));
  const double offset = 2./static_cast<double>(num);
  for (int pt = 0; pt < num; ++pt) {
    const double y = pt*offset - 1 + 0.5*offset;
    const double r = std::sqrt(1. - y*y);
    const double phi = pt*increment;
    (*points)[pt].set_coord(0, radius_*std::cos(phi)*r);
    (*points)[pt].set_coord(1, radius_*y);
    (*points)[pt].set_coord(2, radius_*std::sin(phi)*r);
  }
}

}  // namespace feasst

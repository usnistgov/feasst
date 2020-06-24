#include <cmath>
#include "utils/include/serialize.h"
#include "shape/include/cylinder.h"

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
    const Position point1) {
  class_name_ = "Cylinder";
  args_.init(args);
  radius_ = args_.key("radius").dble();
  point0_ = point0;
  point1_ = point1;
}

double Cylinder::nearest_distance(const Position& point) const {
  return point.nearest_distance_to_axis(point0_, point1_) - radius_;
}

void Cylinder::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_shape_(ostr);
  feasst_serialize_version(629, ostr);
  feasst_serialize(radius_, ostr);
  feasst_serialize_fstobj(point0_, ostr);
  feasst_serialize_fstobj(point1_, ostr);
}

Cylinder::Cylinder(std::istream& istr) : Shape(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(629 == version, version);
  feasst_deserialize(&radius_, istr);
  feasst_deserialize_fstobj(&point0_, istr);
  feasst_deserialize_fstobj(&point1_, istr);
}

}  // namespace feasst

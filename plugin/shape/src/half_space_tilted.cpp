#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "shape/include/half_space_tilted.h"

namespace feasst {

class MapHalfSpaceTilted {
 public:
  MapHalfSpaceTilted() {
    auto obj = MakeHalfSpaceTilted(Position({0, 0, 1}), 1.);
    obj->deserialize_map()["HalfSpaceTilted"] = obj;
  }
};

static MapHalfSpaceTilted mapper_ = MapHalfSpaceTilted();

HalfSpaceTilted::HalfSpaceTilted(const Position& point0,
    const Position& point1) {
  class_name_ = "HalfSpaceTilted";
  ASSERT(std::abs(point0.distance(point1)) > NEAR_ZERO, "point0: "
    << point0.str() << " and point1: " << point1.str() << " are too close.");
  Position unit_normal = point1;
  unit_normal.subtract(point0);
  unit_normal.normalize();
  DEBUG("norm " << unit_normal.str());
  double dist_from_origin;
  const double point1_dist = point1.distance();
  if (std::abs(point1_dist) > NEAR_ZERO) {
    dist_from_origin = point1.dot_product(point0)/point1.distance();
  } else {
    dist_from_origin = -point0.distance();
  }
  DEBUG("dist " << dist_from_origin);
  init_(unit_normal, dist_from_origin);
}

HalfSpaceTilted::HalfSpaceTilted(const Position& unit_normal,
    const double distance_from_origin) {
  init_(unit_normal, distance_from_origin);
}

void HalfSpaceTilted::init_(const Position& unit_normal,
    const double distance_from_origin) {
  class_name_ = "HalfSpaceTilted";
  unit_normal_ = unit_normal;
  distance_from_origin_ = distance_from_origin;
}

// See https://mathworld.wolfram.com/Plane.html
// and https://mathworld.wolfram.com/Point-PlaneDistance.html
// Eq 2, for a plane with normal vector n=(a,b,c) through point x0
// n.(x-x0) = 0 leads to ax+by+cz+d=0 where d=-ax0-by0-cz0=-n.x0
// The distance of the plane from the origin, p=d/|n|
// The signed distance D=|n|.x0+p is positive if x0 is on side of n
double HalfSpaceTilted::nearest_distance(const Position& point) const {
  DEBUG("point " << point.str());
  DEBUG("norm " << unit_normal_.str());
  const double dist = -unit_normal_.dot_product(point) + distance_from_origin_;
  DEBUG(unit_normal_.dot_product(point));
  DEBUG("dist " << dist);
  return dist;
}

void HalfSpaceTilted::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_half_space_tilted_(ostr);
}

void HalfSpaceTilted::serialize_half_space_tilted_(std::ostream& ostr) const {
  serialize_shape_(ostr);
  feasst_serialize_version(7690, ostr);
  feasst_serialize(distance_from_origin_, ostr);
  feasst_serialize_fstobj(unit_normal_, ostr);
}

HalfSpaceTilted::HalfSpaceTilted(std::istream& istr) : Shape(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(7690 == version, version);
  feasst_deserialize(&distance_from_origin_, istr);
  feasst_deserialize_fstobj(&unit_normal_, istr);
}

}  // namespace feasst

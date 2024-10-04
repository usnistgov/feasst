#include <cmath>
#include "utils/include/arguments_extra.h"
#include "utils/include/serialize.h"
#include "shape/include/cylinder.h"

namespace feasst {

FEASST_MAPPER(Cylinder, argtype({{"radius", "1"}, {"first_point", "f"},
  {"f0", "0"}, {"second_point", "s"}, {"s0", "0"}}));

Cylinder::Cylinder(argtype * args) {
  class_name_ = "Cylinder";
  radius_ = dble("radius", args);
  point0_ = Position(parse_dimensional(str("first_point", args), args, 4));
  point1_ = Position(parse_dimensional(str("second_point", args), args, 4));
}
Cylinder::Cylinder(argtype args) : Cylinder(&args) {
  feasst_check_all_used(args);
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

#include <cmath>
#include "utils/include/arguments_extra.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/supertoroid.h"

namespace feasst {

class MapSupertoroid {
 public:
  MapSupertoroid() {
    auto obj = MakeSupertoroid();
    obj->deserialize_map()["Supertoroid"] = obj;
  }
};

static MapSupertoroid mapper_ = MapSupertoroid();

Supertoroid::Supertoroid(argtype * args) : Shape() {
  class_name_ = "Supertoroid";
  a1_ = dble("a1", args, 1);
  a2_ = dble("a2", args, 1);
  a3_ = dble("a3", args, 1);
  a4_ = dble("a4", args, 0);
  epsilon1_ = dble("epsilon1", args, 1);
  epsilon2_ = dble("epsilon2", args, 1);
  if (used("center", *args)) {
    center_ = Position(parse_dimensional(str("center", args), args, 4));
  } else {
    center_ = Position({0, 0, 0});
  }
}
Supertoroid::Supertoroid(argtype args) : Supertoroid(&args) {
  feasst_check_all_used(args);
}

double Supertoroid::nearest_distance(const Position& point) const {
  FATAL("not implemented");
}

bool Supertoroid::is_inside(const Position& point) const {
  double f = std::pow(point.coord(0)/a1_, 2./epsilon2_);
  f += std::pow(point.coord(1)/a2_, 2./epsilon2_);
  f = std::pow(f, epsilon2_/2.);
  f = std::pow(f - a4_, 2./epsilon1_);
  f += std::pow(point.coord(2)/a3_, 2./epsilon1_);
  if (f < 1.) {
    return true;
  }
  return false;
}

bool Supertoroid::is_inside(const Position& point, const double diameter) const {
  return is_inside(point);
}

void Supertoroid::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_shape_(ostr);
  feasst_serialize_version(5635, ostr);
  feasst_serialize(a1_, ostr);
  feasst_serialize(a2_, ostr);
  feasst_serialize(a3_, ostr);
  feasst_serialize(a4_, ostr);
  feasst_serialize(epsilon1_, ostr);
  feasst_serialize(epsilon2_, ostr);
  feasst_serialize_fstobj(center_, ostr);
}

Supertoroid::Supertoroid(std::istream& istr) : Shape(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(5635 == version, version);
  feasst_deserialize(&a1_, istr);
  feasst_deserialize(&a2_, istr);
  feasst_deserialize(&a3_, istr);
  feasst_deserialize(&a4_, istr);
  feasst_deserialize(&epsilon1_, istr);
  feasst_deserialize(&epsilon2_, istr);
  feasst_deserialize_fstobj(&center_, istr);
}

double Supertoroid::surface_area() const { FATAL("not implemented"); }

double Supertoroid::volume() const { FATAL("not implemented"); }

}  // namespace feasst

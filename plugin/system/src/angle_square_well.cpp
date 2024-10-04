#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "configuration/include/bond.h"
#include "system/include/angle_square_well.h"

namespace feasst {

FEASST_MAPPER(AngleSquareWell,);

std::shared_ptr<BondThreeBody> AngleSquareWell::create(std::istream& istr) const {
  return std::make_shared<AngleSquareWell>(istr);
}

AngleSquareWell::AngleSquareWell(std::istream& istr) : BondThreeBody(istr) {
  // ASSERT(class_name_ == "AngleSquareWell", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(846 == version, "mismatch version: " << version);
}

void AngleSquareWell::serialize_angle_square_well_(std::ostream& ostr) const {
  serialize_bond_three_body_(ostr);
  feasst_serialize_version(846, ostr);
}

void AngleSquareWell::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_angle_square_well_(ostr);
}

double AngleSquareWell::energy(const double radians, const Bond& angle) const {
  const double minimum = degrees_to_radians(angle.property("minimum"));
  const double maximum = degrees_to_radians(angle.property("maximum"));
  TRACE("radians " << radians);
  ASSERT(!std::isnan(radians), "radians is nan");
  if (radians < minimum || radians > maximum) {
    return NEAR_INFINITY;
  }
  return 0.;
}

}  // namespace feasst

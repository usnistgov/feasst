#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "system/include/angle_square_well.h"

namespace feasst {

class MapAngleSquareWell {
 public:
  MapAngleSquareWell() {
    auto obj = MakeAngleSquareWell();
    obj->deserialize_map()["AngleSquareWell"] = obj;
  }
};

static MapAngleSquareWell mapper_ = MapAngleSquareWell();

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

double AngleSquareWell::random_angle_radians(const Angle& angle, const double beta,
    const int dimension, Random * random) const {
  return degrees_to_radians(random->uniform_real(angle.property("minimum"),
                                                 angle.property("maximum")));
}

}  // namespace feasst

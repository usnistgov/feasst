#include "utils/include/serialize.h"
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

}  // namespace feasst

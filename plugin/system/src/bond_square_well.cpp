#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "system/include/bond_square_well.h"

namespace feasst {

class MapBondSquareWell {
 public:
  MapBondSquareWell() {
    auto obj = MakeBondSquareWell();
    obj->deserialize_map()["BondSquareWell"] = obj;
  }
};

static MapBondSquareWell mapper_ = MapBondSquareWell();

std::shared_ptr<BondTwoBody> BondSquareWell::create(std::istream& istr) const {
  return std::make_shared<BondSquareWell>(istr);
}

BondSquareWell::BondSquareWell(std::istream& istr) : BondTwoBody(istr) {
  // ASSERT(class_name_ == "BondSquareWell", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(344 == version, "mismatch version: " << version);
}

void BondSquareWell::serialize_bond_square_well_(std::ostream& ostr) const {
  serialize_bond_two_body_(ostr);
  feasst_serialize_version(344, ostr);
}

void BondSquareWell::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bond_square_well_(ostr);
}

double BondSquareWell::energy(
    const Position& relative,
    const Bond& bond) const {
  const double length = bond.property("length");
  const double delta = bond.property("delta");
  if (std::abs(relative.distance() - length) > 0.5*delta) {
    return NEAR_INFINITY;
  }
  return 0.;
}

}  // namespace feasst

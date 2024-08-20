#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/bond.h"
#include "system/include/rigid_bond.h"

namespace feasst {

class MapRigidBond {
 public:
  MapRigidBond() {
    auto obj = MakeRigidBond();
    obj->deserialize_map()["RigidBond"] = obj;
  }
};

static MapRigidBond mapper_ = MapRigidBond();

std::shared_ptr<BondTwoBody> RigidBond::create(std::istream& istr) const {
  return std::make_shared<RigidBond>(istr);
}

RigidBond::RigidBond(std::istream& istr) : BondTwoBody(istr) {
  // ASSERT(class_name_ == "RigidBond", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7690 == version, "mismatch version: " << version);
}

void RigidBond::serialize_rigid_bond_(std::ostream& ostr) const {
  serialize_bond_two_body_(ostr);
  feasst_serialize_version(7690, ostr);
}

void RigidBond::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_rigid_bond_(ostr);
}

double RigidBond::energy(const double distance, const Bond& bond) const {
  const double length = bond.property("length");
  const double delta = bond.property("delta");
  if (std::abs(distance - length) > delta) {
    FATAL("distance(" << distance << ")-length(" << length << ")=" <<
      distance - length << " > delta: " << delta);
  }
  return 0.;
}

double RigidBond::random_distance(const Bond& bond, const double beta,
    const int dimen, Random * random) const {
  return bond.property("length");
}

}  // namespace feasst

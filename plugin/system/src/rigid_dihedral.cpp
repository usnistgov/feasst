#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "configuration/include/bond.h"
#include "system/include/rigid_dihedral.h"

namespace feasst {

FEASST_MAPPER(RigidDihedral,);

std::shared_ptr<BondFourBody> RigidDihedral::create(std::istream& istr) const {
  return std::make_shared<RigidDihedral>(istr);
}

RigidDihedral::RigidDihedral(std::istream& istr) : BondFourBody(istr) {
  // ASSERT(class_name_ == "RigidDihedral", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(1274 == version, "mismatch version: " << version);
}

void RigidDihedral::serialize_rigid_dihedral_(std::ostream& ostr) const {
  serialize_bond_four_body_(ostr);
  feasst_serialize_version(1274, ostr);
}

void RigidDihedral::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_rigid_dihedral_(ostr);
}

double RigidDihedral::energy(const double radians, const Bond& dihedral) const {
  TRACE("radians " << radians);
  const double theta = degrees_to_radians(dihedral.property("degrees"));
  const double delta = degrees_to_radians(dihedral.property("delta"));
  ASSERT(!std::isnan(radians), "radians is nan");
  if (std::abs(radians - theta) > delta) {
    return NEAR_INFINITY;
  }
  return 0.;
}

double RigidDihedral::random_dihedral_radians(const Dihedral& dihedral,
    const double beta, const int dimension, Random * random) const {
  return degrees_to_radians(dihedral.property("degrees"));
}

}  // namespace feasst

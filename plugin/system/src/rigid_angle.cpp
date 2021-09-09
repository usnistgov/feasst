#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "system/include/rigid_angle.h"

namespace feasst {

class MapRigidAngle {
 public:
  MapRigidAngle() {
    auto obj = MakeRigidAngle();
    obj->deserialize_map()["RigidAngle"] = obj;
  }
};

static MapRigidAngle mapper_ = MapRigidAngle();

std::shared_ptr<BondThreeBody> RigidAngle::create(std::istream& istr) const {
  return std::make_shared<RigidAngle>(istr);
}

RigidAngle::RigidAngle(std::istream& istr) : BondThreeBody(istr) {
  // ASSERT(class_name_ == "RigidAngle", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(5968 == version, "mismatch version: " << version);
}

void RigidAngle::serialize_rigid_angle_(std::ostream& ostr) const {
  serialize_bond_three_body_(ostr);
  feasst_serialize_version(5968, ostr);
}

void RigidAngle::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_rigid_angle_(ostr);
}

double RigidAngle::energy(const double radians, const Bond& angle) const {
  const double theta = degrees_to_radians(angle.property("degrees"));
  const double delta = degrees_to_radians(angle.property("delta"));
  ASSERT(!std::isnan(radians), "radians is nan");
  if (std::abs(radians - theta) > delta) {
    FATAL("radians(" << radians << ")-theta(" << theta << ")=" << radians - theta
      << " > delta: " << delta);
  }
  return 0.;
}

double RigidAngle::random_angle_radians(const Angle& angle, const double beta,
    const int dimension, Random * random) const {
  return degrees_to_radians(angle.property("degrees"));
}

void RigidAngle::random_branch(
    const Angle& a2a1m1,
    const Angle& a2a1m2,
    const Angle& m1a1m2,
    const double beta,
    double * radians_a2a1m1,
    double * radians_a2a1m2,
    double * radians_m1a1m2,
    Random * random) const {
  ASSERT(a2a1m1.model() == "RigidAngle" &&
         a2a1m2.model() == "RigidAngle" &&
         m1a1m2.model() == "RigidAngle", "Branch model mismatch");
  *radians_a2a1m1 = random_angle_radians(a2a1m1, beta, 3, random);
  *radians_a2a1m2 = random_angle_radians(a2a1m2, beta, 3, random);
  *radians_m1a1m2 = random_angle_radians(m1a1m2, beta, 3, random);
}

}  // namespace feasst

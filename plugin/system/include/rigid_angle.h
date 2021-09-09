
#ifndef FEASST_SYSTEM_RIGID_ANGLE_H_
#define FEASST_SYSTEM_RIGID_ANGLE_H_

#include <memory>
#include "system/include/bond_three_body.h"

namespace feasst {

/**
  U(r) = 0 when the absolute value of the difference between the actual angle
  and the degrees parameter is less than delta.
  Otherwise, U(r) = NEAR_INFINITY.
 */
class RigidAngle : public BondThreeBody {
 public:
  RigidAngle() {}
  double energy(const double radians, const Bond& angle) const override;
  double random_angle_radians(const Angle& angle, const double beta,
    const int dimension, Random * random) const override;
  void random_branch(
    const Angle& a2a1m1,
    const Angle& a2a1m2,
    const Angle& m1a1m2,
    const double beta,
    double * radians_a2a1m1,
    double * radians_a2a1m2,
    double * radians_m1a1m2,
    Random * random) const override;
  std::shared_ptr<BondThreeBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit RigidAngle(std::istream& istr);
  virtual ~RigidAngle() {}

 protected:
  void serialize_rigid_angle_(std::ostream& ostr) const;
};

inline std::shared_ptr<RigidAngle> MakeRigidAngle() {
  return std::make_shared<RigidAngle>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_RIGID_ANGLE_H_


#ifndef FEASST_SYSTEM_RIGID_BOND_H_
#define FEASST_SYSTEM_RIGID_BOND_H_

#include <memory>
#include "system/include/bond_two_body.h"

namespace feasst {

/**
  U(r) = 0 when the absolute value of the difference between the actual length
  and the length parameter is less than delta.
  Otherwise, U(r) = NEAR_INFINITY.
 */
class RigidBond : public BondTwoBody {
 public:
  RigidBond() {}
  double energy(const double distance, const Bond& bond) const override;
  double random_distance(const Bond& bond, const double beta, const int dimen,
    Random * random) const override;
  std::shared_ptr<BondTwoBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit RigidBond(std::istream& istr);
  virtual ~RigidBond() {}

 protected:
  void serialize_rigid_bond_(std::ostream& ostr) const;
};

inline std::shared_ptr<RigidBond> MakeRigidBond() {
  return std::make_shared<RigidBond>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_RIGID_BOND_H_

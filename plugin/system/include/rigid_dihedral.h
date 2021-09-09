
#ifndef FEASST_SYSTEM_RIGID_DIHEDRAL_H_
#define FEASST_SYSTEM_RIGID_DIHEDRAL_H_

#include <memory>
#include "system/include/bond_four_body.h"

namespace feasst {

/**
  U(r) = 0 when the absolute value of the difference between the actual dihedral
  and the degrees parameter is less than delta.
  Otherwise, U(r) = NEAR_INFINITY.
 */
class RigidDihedral : public BondFourBody {
 public:
  RigidDihedral() {}
  double energy(const double radians, const Bond& dihedral) const override;
  double random_dihedral_radians(const Dihedral& dihedral, const double beta,
    const int dimension, Random * random) const override;
  std::shared_ptr<BondFourBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit RigidDihedral(std::istream& istr);
  virtual ~RigidDihedral() {}

 protected:
  void serialize_rigid_dihedral_(std::ostream& ostr) const;
};

inline std::shared_ptr<RigidDihedral> MakeRigidDihedral() {
  return std::make_shared<RigidDihedral>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_RIGID_DIHEDRAL_H_

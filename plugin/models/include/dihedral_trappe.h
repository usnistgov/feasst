
#ifndef FEASST_MODELS_DIHEDRAL_TRAPPE_H_
#define FEASST_MODELS_DIHEDRAL_TRAPPE_H_

#include <memory>
#include "system/include/bond_four_body.h"

namespace feasst {

/**
  See the Dihedral section of http://trappe.oit.umn.edu/

  \f$U(\phi)=c_0+c_1[1+\cos(\phi)]+c_2[1-\cos(2\phi)]+c_3[1+\cos(3\phi)]\f$
 */
class DihedralTraPPE : public BondFourBody {
 public:
  explicit DihedralTraPPE(const argtype& args = argtype()) {}
  double energy(const double radians, const Bond& dihedral) const override;
  std::shared_ptr<BondFourBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit DihedralTraPPE(std::istream& istr);
  virtual ~DihedralTraPPE() {}

 protected:
  void serialize_dihedral_trappe_(std::ostream& ostr) const;
};

inline std::shared_ptr<DihedralTraPPE> MakeDihedralTraPPE(
    const argtype &args = argtype()) {
  return std::make_shared<DihedralTraPPE>(args);
}

}  // namespace feasst

#endif  // FEASST_MODELS_DIHEDRAL_TRAPPE_H_

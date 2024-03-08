
#ifndef FEASST_MODELS_DIHEDRAL_RYCKAERT_BELLEMANS_H_
#define FEASST_MODELS_DIHEDRAL_RYCKAERT_BELLEMANS_H_

#include <memory>
#include "system/include/bond_four_body.h"

namespace feasst {

/**
  See the following manuscripts:
  https://doi.org/10.1039/DC9786600095

  As described in the L-OPLS implementation:
  https://doi.org/10.1021/ct200908r

  \f$U(\phi)=c_0+c_1(\cos\phi)+c_2(\cos\phi)^2+c_3(\cos\phi)^3\f$
 */
class DihedralRyckaertBellemans : public BondFourBody {
 public:
  explicit DihedralRyckaertBellemans(const argtype& args = argtype()) {}
  double energy(const double radians, const Bond& dihedral) const override;
  std::shared_ptr<BondFourBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit DihedralRyckaertBellemans(std::istream& istr);
  virtual ~DihedralRyckaertBellemans() {}

 protected:
  void serialize_dihedral_ryckaert_bellemans_(std::ostream& ostr) const;
};

inline std::shared_ptr<DihedralRyckaertBellemans> MakeDihedralRyckaertBellemans(
    const argtype &args = argtype()) {
  return std::make_shared<DihedralRyckaertBellemans>(args);
}

}  // namespace feasst

#endif  // FEASST_MODELS_DIHEDRAL_RYCKAERT_BELLEMANS_H_

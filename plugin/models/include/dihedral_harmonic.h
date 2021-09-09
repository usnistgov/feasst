
#ifndef FEASST_MODELS_DIHEDRAL_HARMONIC_H_
#define FEASST_MODELS_DIHEDRAL_HARMONIC_H_

#include <memory>
#include "system/include/bond_four_body.h"

namespace feasst {

/**
  U(angle) = k_energy_per_radian_sq*(angle - equilibrium_degrees)^2
  with parameters given in Dihedral Properties.

  The usual 1/2 factor is not included, but can be incorporated into
  the k parameter manually by the user input to the forcefield file.
 */
class DihedralHarmonic : public BondFourBody {
 public:
  DihedralHarmonic() {}
  double energy(const double radians, const Bond& dihedral) const override;
  std::shared_ptr<BondFourBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit DihedralHarmonic(std::istream& istr);
  virtual ~DihedralHarmonic() {}

 protected:
  void serialize_dihedral_harmonic_(std::ostream& ostr) const;
};

inline std::shared_ptr<DihedralHarmonic> MakeDihedralHarmonic() {
  return std::make_shared<DihedralHarmonic>();
}

}  // namespace feasst

#endif  // FEASST_MODELS_DIHEDRAL_HARMONIC_H_

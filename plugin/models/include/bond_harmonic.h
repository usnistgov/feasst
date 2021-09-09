
#ifndef FEASST_MODELS_BOND_HARMONIC_H_
#define FEASST_MODELS_BOND_HARMONIC_H_

#include <memory>
#include "system/include/bond_two_body.h"

namespace feasst {

/**
  U(r) = k_energy_per_length_sq*(r - equilibrium_length)^2
  with parameters given in Bond Properties.

  The usual 1/2 factor is not included, but can be incorporated into
  the k parameter manually by the user input to the forcefield file.
 */
class BondHarmonic : public BondTwoBody {
 public:
  BondHarmonic() {}
  double energy(const double distance, const Bond& bond) const override;

  /**
    Return a randomly selected bond length with harmonic potential as described
    in Frenkel and Smit, Alg 43, page 578 and Allen and Tildesley, Section G.3.
   
    \f$P(l)dl \propto l**2\exp[-\beta U(length)]dl\f$

    The maximal length is 3 sigma beyond the mean.
    If 2D, use the accept-reject method
   */
  double random_distance(const Bond& bond, const double beta, const int dimen,
    Random * random) const override;
  std::shared_ptr<BondTwoBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit BondHarmonic(std::istream& istr);
  virtual ~BondHarmonic() {}

 protected:
  void serialize_bond_harmonic_(std::ostream& ostr) const;
};

inline std::shared_ptr<BondHarmonic> MakeBondHarmonic() {
  return std::make_shared<BondHarmonic>();
}

}  // namespace feasst

#endif  // FEASST_MODELS_BOND_HARMONIC_H_

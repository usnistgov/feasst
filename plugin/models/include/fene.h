
#ifndef FEASST_MODELS_FENE_H_
#define FEASST_MODELS_FENE_H_

#include <memory>
#include "system/include/bond_two_body.h"

namespace feasst {

/**
  FENE is a two-body finite extensible nonlinear elastic two-body potential.
  R. B. Bird et al., Dynamics of Polymeric Liquids (Wiley, New York, 1977), Vols. 1 and 2.

  \f$ U = -0.5*\f$k_energy_per_length_sq\f$*\f$R_0\f$^2\ln\left[1-\left(\frac{r}{R_0}\right)^2\right] \f$

  where k_energy_per_length_sq and R_0 are input bond parameters in the forcefield.

  Note that FENE is typically coupled with a WCA potential as described in
  LennardJonesCutShift.
 */
class FENE : public BondTwoBody {
 public:
  FENE() {}
  double energy_fene(const double distance, const double k, const double R0) const;
  double energy(const double distance, const Bond& bond) const override;
  double random_distance(const Bond& bond, const double beta,
    const int dimen, Random * random) const override;
  std::shared_ptr<BondTwoBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit FENE(std::istream& istr);
  virtual ~FENE() {}

 protected:
  void serialize_fene_(std::ostream& ostr) const;
};

inline std::shared_ptr<FENE> MakeFENE() { return std::make_shared<FENE>(); }

}  // namespace feasst

#endif  // FEASST_MODELS_FENE_H_

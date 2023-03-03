
#ifndef FEASST_CHARGE_DEBYE_HUCKEL_H_
#define FEASST_CHARGE_DEBYE_HUCKEL_H_

#include "system/include/model_two_body.h"

namespace feasst {

/**
  Compute energy between two point charges, \f$q_i\f$ and \f$q_j\f$ with a
  Debye-Huckel potential.

  \f$U = q_i q_j \chi \exp(-\kappa r) / \epsilon r\f$

  where \f$r\f$ is the separation distance,
  \f$\kappa\f$ is the inverse of the Debye screening length,
  \f$\epsilon\f$ is the dilectric constant (dimensionless),
  and \f$\chi\f$ is the charge conversion factor assuming the following units:
    1. length: Angstroms
    2. energy: kJ/mol
    3. charge: elementary

  An optional smoothing_distance linearly interpolates the energy to zero at
  the cutoff starting at a distance of cutoff - smoothing_distance.

  Avoid singularity by returning a large, positive number when \f$r\f$ is near zero.
 */
class DebyeHuckel : public ModelTwoBody {
 public:
  /**
    args:
    - kappa: as described above.
    - dielectric: as described above.
    - smoothing_distance: as described above.
      Disabled when negative (default: -1).
   */
  explicit DebyeHuckel(argtype args = argtype());
  explicit DebyeHuckel(argtype * args);

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  void precompute(const ModelParams& existing) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<DebyeHuckel>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<DebyeHuckel>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit DebyeHuckel(std::istream& istr);
  virtual ~DebyeHuckel() {}

 private:
  double conversion_factor_ = 0.;
  double kappa_;
  double dielectric_;
  double smoothing_distance_;
};

inline std::shared_ptr<DebyeHuckel> MakeDebyeHuckel(
    argtype args = argtype()) {
  return std::make_shared<DebyeHuckel>(args);
}

}  // namespace feasst

#endif  // FEASST_CHARGE_DEBYE_HUCKEL_H_

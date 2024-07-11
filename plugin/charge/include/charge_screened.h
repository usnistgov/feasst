
#ifndef FEASST_CHARGE_CHARGE_SCREENED_H_
#define FEASST_CHARGE_CHARGE_SCREENED_H_

#include "system/include/model_two_body.h"

namespace feasst {

class Table1D;

/**
  Compute energy between two point charges, \f$q_i\f$ and \f$q_j\f$ with a
  Gaussian screening cloud as utilized by the Ewald summation.

  \f$U = q_i q_j \chi erfc(\alpha r)/r\f$

  where \f$erfc\f$ is the complimentary error function,
  \f$r\f$ is the separation distance,
  and \f$\chi\f$ is the charge conversion factor assuming the following units:
    1. length: Angstroms
    2. energy: kJ/mol
    3. charge: elementary

  Avoid Coulomb explosion by returning a large number when \f$r\f$ is near zero.
 */
class ChargeScreened : public ModelTwoBody {
 public:
  //@{
  /** @name Arguments
    - hard_sphere_threshold: return NEAR_INFINITY when distance is less than
      this threshold (default: 0.2).
    - erfc_table_size: size of forward-difference-interpolated tabular potential
      (default: 0).
      Disable table if this value is less than or equal to zero.
      This table does not depend on site type, so may use less memory than the
      generate Potential::table_size.
   */
  explicit ChargeScreened(argtype args = argtype());
  explicit ChargeScreened(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  void precompute(const ModelParams& existing) override;

  /// Return the erfc table.
  const Table1D * erfc(const double distance_squared) const {
    return erfc_.get(); }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ChargeScreened>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ChargeScreened>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ChargeScreened(std::istream& istr);
  virtual ~ChargeScreened() {}

  //@}
 private:
  double alpha_;
  double conversion_factor_;
  double hard_sphere_threshold_sq_;
  int erfc_table_size_;
  std::shared_ptr<Table1D> erfc_;
  void init_erfc_(const double cutoff);
};

inline std::shared_ptr<ChargeScreened> MakeChargeScreened(
    argtype args = argtype()) {
  return std::make_shared<ChargeScreened>(args);
}

}  // namespace feasst

#endif  // FEASST_CHARGE_CHARGE_SCREENED_H_

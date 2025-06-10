
#ifndef FEASST_CHARGE_CHARGE_SCREENED_INTRA_H_
#define FEASST_CHARGE_CHARGE_SCREENED_INTRA_H_

#include "system/include/model_two_body.h"

namespace feasst {

/**
  Compute the correction energy used to remove the spurious intra-particle
  interactions that the Ewald summation includes

  \f$U = -q_i q_j \chi erf(\alpha r)/r\f$

  see ChargeScreened for details.
 */
class ChargeScreenedIntra : public ModelTwoBody {
 public:
  ChargeScreenedIntra() { class_name_ = "ChargeScreenedIntra"; }

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  void precompute(const Configuration& config) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ChargeScreenedIntra>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ChargeScreenedIntra>(); }
  void serialize(std::ostream& ostr) const override;
  explicit ChargeScreenedIntra(std::istream& istr);
  virtual ~ChargeScreenedIntra() {}

 private:
  double alpha_;
  double conversion_factor_;
};

inline std::shared_ptr<ChargeScreenedIntra> MakeChargeScreenedIntra() {
  return std::make_shared<ChargeScreenedIntra>();
}

}  // namespace feasst

#endif  // FEASST_CHARGE_CHARGE_SCREENED_INTRA_H_

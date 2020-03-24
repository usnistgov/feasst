
#ifndef FEASST_EWALD_CHARGE_SCREENED_INTRA_H_
#define FEASST_EWALD_CHARGE_SCREENED_INTRA_H_

#include "system/include/model_two_body.h"

namespace feasst {

/**
  Assumes units:
    1. length: Angstroms
    2. energy: kJ/mol
    3. charge: elementary
 */
class ChargeScreenedIntra : public ModelTwoBody {
 public:
  ChargeScreenedIntra() {}

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override;

  void precompute(const ModelParams& existing) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ChargeScreenedIntra>(istr); }

  void serialize(std::ostream& ostr) const override;
  explicit ChargeScreenedIntra(std::istream& istr);
  virtual ~ChargeScreenedIntra() {}

 private:
  const std::string class_name_ = "ChargeScreenedIntra";
  double alpha_;
  double conversion_factor_;
};

inline std::shared_ptr<ChargeScreenedIntra> MakeChargeScreenedIntra() {
  return std::make_shared<ChargeScreenedIntra>();
}

}  // namespace feasst

#endif  // FEASST_EWALD_CHARGE_SCREENED_INTRA_H_

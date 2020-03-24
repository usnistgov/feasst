
#ifndef FEASST_EWALD_CHARGE_SCREENED_H_
#define FEASST_EWALD_CHARGE_SCREENED_H_

#include "system/include/model_two_body.h"

namespace feasst {

/**
  Assumes units:
    1. length: Angstroms
    2. energy: kJ/mol
    3. charge: elementary
 */
class ChargeScreened : public ModelTwoBody {
 public:
  ChargeScreened() {}

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override;

  void precompute(const ModelParams& existing) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ChargeScreened>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit ChargeScreened(std::istream& istr);
  virtual ~ChargeScreened() {}

 private:
  const std::string class_name_ = "ChargeScreened";
  double alpha_;
  double conversion_factor_;
};

inline std::shared_ptr<ChargeScreened> MakeChargeScreened() {
  return std::make_shared<ChargeScreened>();
}

}  // namespace feasst

#endif  // FEASST_EWALD_CHARGE_SCREENED_H_

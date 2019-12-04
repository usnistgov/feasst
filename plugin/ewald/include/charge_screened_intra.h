
#ifndef FEASST_EWALD_CHARGE_SCREENED_INTRA_H_
#define FEASST_EWALD_CHARGE_SCREENED_INTRA_H_

#include "system/include/model_two_body.h"
#include "math/include/constants.h"
#include "system/include/physical_constants.h"

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
      const ModelParams& model_params) const override {
    const double mixed_charge = model_params.mixed_charge()[type1][type2];
    const double distance = std::sqrt(squared_distance);
    return -mixed_charge*charge_conversion*erf(alpha_*distance)/distance;
  }

  void precompute(const ModelParams& existing) override {
    alpha_ = existing.property("alpha");
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ChargeScreenedIntra>(istr); }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
    feasst_serialize(alpha_, ostr);
  }

  ChargeScreenedIntra(std::istream& istr) {
    feasst_deserialize_version(istr);
    feasst_deserialize(&alpha_, istr);
  }

  virtual ~ChargeScreenedIntra() {}

 private:
  const std::string class_name_ = "ChargeScreenedIntra";
  double alpha_;
};

inline std::shared_ptr<ChargeScreenedIntra> MakeChargeScreenedIntra() {
  return std::make_shared<ChargeScreenedIntra>();
}

}  // namespace feasst

#endif  // FEASST_EWALD_CHARGE_SCREENED_INTRA_H_

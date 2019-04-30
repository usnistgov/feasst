
#ifndef FEASST_EWALD_MODEL_CHARGE_SCREENED_H_
#define FEASST_EWALD_MODEL_CHARGE_SCREENED_H_

#include "core/include/model_two_body.h"
#include "core/include/constants.h"
#include "core/include/physical_constants.h"

namespace feasst {

/**
  Assumes units:
    1. length: Angstroms
    2. energy: kJ/mol
    3. charge: elementary
 */
class ModelChargeScreened : public ModelTwoBody {
 public:
  ModelChargeScreened() {}

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override {
    const double mixed_charge = model_params.mixed_charge()[type1][type2];
    const double distance = std::sqrt(squared_distance);
    return mixed_charge*charge_conversion*erfc(alpha_*distance)/distance;
  }

  void precompute(const ModelParams& existing) override {
    alpha_ = existing.property("alpha");
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelChargeScreened>(istr);
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize(alpha_, ostr);
    feasst_serialize_version(1, ostr);
  }

  ModelChargeScreened(std::istream& istr) {
    feasst_deserialize_version(istr);
    feasst_deserialize(&alpha_, istr);
  }

  virtual ~ModelChargeScreened() {}

 private:
  const std::string class_name_ = "ModelChargeScreened";
  double alpha_;
};

inline std::shared_ptr<ModelChargeScreened> MakeModelChargeScreened() {
  return std::make_shared<ModelChargeScreened>();
}

}  // namespace feasst

#endif  // FEASST_EWALD_MODEL_CHARGE_SCREENED_H_

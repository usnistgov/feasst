
#ifndef FEASST_EWALD_MODEL_CHARGE_SELF_H_
#define FEASST_EWALD_MODEL_CHARGE_SELF_H_

#include "system/include/model_one_body.h"
#include "math/include/constants.h"
#include "system/include/physical_constants.h"

namespace feasst {

/**
  Assumes units:
    1. length: Angstroms
    2. energy: kJ/mol
    3. charge: elementary
 */
class ModelChargeSelf : public ModelOneBody {
 public:
  ModelChargeSelf() {}

  double energy(
      const Site& site,
      const Configuration * config,
      const ModelParams& model_params) const override {
    const int type = site.type();
    const double charge = model_params.charge().value(type);
    return -charge*charge*charge_conversion*alpha_/std::sqrt(PI);
  }

  void precompute(const ModelParams& existing) override {
    alpha_ = existing.property("alpha");
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelChargeSelf>(istr); }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
    feasst_serialize(alpha_, ostr);
  }

  ModelChargeSelf(std::istream& istr) {
    feasst_deserialize_version(istr);
    feasst_deserialize(&alpha_, istr);
  }

  virtual ~ModelChargeSelf() {}

 private:
  const std::string class_name_ = "ModelChargeSelf";
  double alpha_;
};

inline std::shared_ptr<ModelChargeSelf> MakeModelChargeSelf() {
  return std::make_shared<ModelChargeSelf>();
}

}  // namespace feasst

#endif  // FEASST_EWALD_MODEL_CHARGE_SELF_H_

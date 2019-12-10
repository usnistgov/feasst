
#ifndef FEASST_EWALD_CHARGE_SELF_H_
#define FEASST_EWALD_CHARGE_SELF_H_

#include "system/include/model_one_body.h"
#include "math/include/constants.h"

namespace feasst {

/**
  Assumes units:
    1. length: Angstroms
    2. energy: kJ/mol
    3. charge: elementary
 */
class ChargeSelf : public ModelOneBody {
 public:
  ChargeSelf() {}

  double energy(
      const Site& site,
      const Configuration * config,
      const ModelParams& model_params) const override {
    const int type = site.type();
    const double charge = model_params.charge().value(type);
    return -charge*charge*conversion_factor_*alpha_/std::sqrt(PI);
  }

  void precompute(const ModelParams& existing) override {
    alpha_ = existing.property("alpha");
    conversion_factor_ = existing.constants()->charge_conversion();
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ChargeSelf>(istr); }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
    feasst_serialize(alpha_, ostr);
    feasst_serialize(conversion_factor_, ostr);
  }

  ChargeSelf(std::istream& istr) {
    feasst_deserialize_version(istr);
    feasst_deserialize(&alpha_, istr);
    feasst_deserialize(&conversion_factor_, istr);
  }

  virtual ~ChargeSelf() {}

 private:
  const std::string class_name_ = "ChargeSelf";
  double alpha_;
  double conversion_factor_;
};

inline std::shared_ptr<ChargeSelf> MakeChargeSelf() {
  return std::make_shared<ChargeSelf>();
}

}  // namespace feasst

#endif  // FEASST_EWALD_CHARGE_SELF_H_


#ifndef FEASST_EWALD_CHARGE_SELF_H_
#define FEASST_EWALD_CHARGE_SELF_H_

#include "system/include/model_one_body.h"

namespace feasst {

/**
  Compute the correction energy used to remove the spurious self-self
  interaction that the Ewald summation includes

  \f$U = -q_i q_i \chi \frac{\alpha}{\sqrt{\pi}} \f$

  see ChargeScreened for details.
 */
class ChargeSelf : public ModelOneBody {
 public:
  ChargeSelf() {}

  double energy(
      const Site& site,
      const Configuration& config,
      const ModelParams& model_params) const override;

  void precompute(const ModelParams& existing) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ChargeSelf>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit ChargeSelf(std::istream& istr);
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

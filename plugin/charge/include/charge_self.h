
#ifndef FEASST_CHARGE_CHARGE_SELF_H_
#define FEASST_CHARGE_CHARGE_SELF_H_

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
  ChargeSelf() { class_name_ = "ChargeSelf"; }

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) override;

  void precompute(const Configuration& config) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ChargeSelf>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ChargeSelf>(); }
  void serialize(std::ostream& ostr) const override;
  explicit ChargeSelf(std::istream& istr);
  virtual ~ChargeSelf() {}

 private:
  double alpha_;
  double conversion_factor_;
};

inline std::shared_ptr<ChargeSelf> MakeChargeSelf() {
  return std::make_shared<ChargeSelf>();
}

}  // namespace feasst

#endif  // FEASST_CHARGE_CHARGE_SELF_H_

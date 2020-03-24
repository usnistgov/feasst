
#ifndef FEASST_CHAIN_TRIAL_SELECT_PERTURBED_H_
#define FEASST_CHAIN_TRIAL_SELECT_PERTURBED_H_

#include "monte_carlo/include/trial_select.h"

namespace feasst {

/// Select the perturbed site. Used for selection of first stage with growth
/// expected ensemble. Particular particle, not random particle.
class TrialSelectPerturbed : public TrialSelect {
 public:
  TrialSelectPerturbed(const argtype& args = argtype()) : TrialSelect(args) {
    class_name_ = "TrialSelectPerturbed";
  }

  bool select(const Select& perturbed, System* system, Random * random) override {
    if (perturbed.num_sites() == 0) return false;
    replace_mobile(perturbed, 0, system->configuration());
    mobile_original_ = mobile_;
    return true;
  }

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectPerturbed(std::istream& istr);
  virtual ~TrialSelectPerturbed() {}
};

inline std::shared_ptr<TrialSelectPerturbed> MakeTrialSelectPerturbed(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectPerturbed>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SELECT_PERTURBED_H_

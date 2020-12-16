
#ifndef FEASST_CHAIN_SELECT_PERTURBED_H_
#define FEASST_CHAIN_SELECT_PERTURBED_H_

#include "monte_carlo/include/trial_select.h"

namespace feasst {

// HWH move to growth expanded ensemble
/// Select the perturbed site. Used for selection of first stage with growth
/// expanded ensemble. Particular particle, not random particle.
class SelectPerturbed : public TrialSelect {
 public:
  SelectPerturbed(const argtype& args = argtype()) : TrialSelect(args) {
    class_name_ = "SelectPerturbed";
  }

  bool select(const Select& perturbed, System* system, Random * random) override {
    if (perturbed.num_sites() == 0) return false;
    replace_mobile(perturbed, 0, system->configuration());
    mobile_original_ = mobile_;
    return true;
  }

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectPerturbed(std::istream& istr);
  virtual ~SelectPerturbed() {}
};

inline std::shared_ptr<SelectPerturbed> MakeSelectPerturbed(
    const argtype &args = argtype()) {
  return std::make_shared<SelectPerturbed>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_SELECT_PERTURBED_H_

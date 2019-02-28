
#ifndef FEASST_CORE_TRIAL_REGROW_H_
#define FEASST_CORE_TRIAL_REGROW_H_

#include "core/include/trial_translate.h"
#include "chain/include/perturb_regrow.h"

namespace feasst {

class TrialRegrow : public TrialMove {
 public:
  TrialRegrow() {
    regrow_ = std::make_shared<PerturbRegrow>();
    set_perturb(regrow_);
  }

  void select(System * system) override {
    regrow_->select_random_end_segment_in_particle(group_index(), system->configuration());
  }

  void move(System * system) override {
    regrow_->regrow(system);
  }

  virtual ~TrialRegrow() {}

 private:
  std::shared_ptr<PerturbRegrow> regrow_;
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_REGROW_H_

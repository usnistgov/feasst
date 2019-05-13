
#ifndef FEASST_CHAIN_TRIAL_CRANKSHAFT_H_
#define FEASST_CHAIN_TRIAL_CRANKSHAFT_H_

#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/perturb_translate.h"
#include "math/include/random.h"
#include "utils/include/utils_io.h"

namespace feasst {

class TrialCrankshaft : public TrialRotate {
 public:
  TrialCrankshaft(
    /**
      max_length : maximum length of selected segment. If -1 (default), then
        randomly select all possible lengths.
     */
    const argtype& args = argtype()) : TrialRotate(args) {
    args_.init(args);
    max_length_ = args_.key("max_length").dflt("-1").integer();
    set_recenter(1);
  }

  void select(System * system) override {
    perturb_rotate_->select_random_segment_in_particle(
      group_index(),
      system->configuration(),
      max_length_);
  }

  void move(Criteria * criteria, System * system) override {
    const Position& pivot = perturb_rotate_->selection().site_positions()[0].front();
    {
      Position axis = perturb_rotate_->selection().site_positions()[0].back();
      axis.subtract(pivot);
      axis.normalize();
      const double max_angle = perturb_rotate_->tunable().value();
      const double angle = random_.uniform_real(-max_angle, max_angle);
      rot_mat_.axis_angle(axis, angle);
    }
    perturb_rotate_->rotate_selection(pivot, rot_mat_, system);
  }

  virtual ~TrialCrankshaft() {}

 private:
  int max_length_;
};

inline std::shared_ptr<TrialCrankshaft> MakeTrialCrankshaft(
    const argtype &args = argtype()) {
  return std::make_shared<TrialCrankshaft>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_CRANKSHAFT_H_

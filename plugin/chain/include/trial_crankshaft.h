
#ifndef FEASST_CHAIN_TRIAL_CRANKSHAFT_H_
#define FEASST_CHAIN_TRIAL_CRANKSHAFT_H_

#include "core/include/trial_rotate.h"
#include "core/include/perturb_translate.h"
#include "core/include/random.h"
#include "core/include/utils_io.h"

namespace feasst {

class TrialCrankshaft : public TrialRotate {
 public:
  TrialCrankshaft(const argtype& args = argtype()) : TrialRotate(args) { set_recenter(1); }

  void select(System * system) override {
    perturb_rotate_->select_random_segment_in_particle(group_index(), system->configuration());
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
};

inline std::shared_ptr<TrialCrankshaft> MakeTrialCrankshaft(
    const argtype &args = argtype()) {
  return std::make_shared<TrialCrankshaft>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_CRANKSHAFT_H_

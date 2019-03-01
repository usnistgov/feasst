
#ifndef FEASST_CHAIN_TRIAL_PIVOT_H_
#define FEASST_CHAIN_TRIAL_PIVOT_H_

#include "core/include/trial_rotate.h"
#include "core/include/utils_io.h"

namespace feasst {

class TrialPivot : public TrialRotate {
 public:
  TrialPivot() : TrialRotate() { set_recenter(1); }

  void select(System * system) override {
    perturb_rotate_->select_random_end_segment_in_particle(group_index(), system->configuration());
  }

  void move(System * system) override {
    // the last site in selection is assumed to be the pivot point based on selection
    const Position& pivot = perturb_rotate_->selection().site_positions()[0].back();
    DEBUG("pivot " << pivot.str());
    random_rotation(pivot);
    perturb_rotate_->rotate_selection(pivot, rot_mat_, system);
  }

  virtual ~TrialPivot() {}
};

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_PIVOT_H_

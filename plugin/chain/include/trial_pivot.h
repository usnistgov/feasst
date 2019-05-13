
#ifndef FEASST_CHAIN_TRIAL_PIVOT_H_
#define FEASST_CHAIN_TRIAL_PIVOT_H_

#include "monte_carlo/include/trial_rotate.h"
#include "utils/include/utils_io.h"

namespace feasst {

class TrialPivot : public TrialRotate {
 public:
    /**
      max_length : maximum length of selected segment. If -1 (default), then
        randomly select all possible lengths.
     */
  TrialPivot(const argtype& args = argtype()) : TrialRotate(args) {
    args_.init(args);
    max_length_ = args_.key("max_length").dflt("-1").integer();
    set_recenter(1);
  }

  void select(System * system) override {
    perturb_rotate_->select_random_end_segment_in_particle(
      group_index(),
      system->configuration(),
      max_length_);
  }

  void move(Criteria * criteria, System * system) override {
    // find the pivot index. If first site is end point then pivot is last.
    // otherwise last site is endpoint and pivot is first.
    int pivot_index = 0;
    if (perturb_rotate_->selection().site_indices()[0][0] == 0) {
      pivot_index = perturb_rotate_->selection().num_sites() - 1;
    }
    const Position& pivot = perturb_rotate_->selection().site_positions()[0][pivot_index];
    DEBUG("pivot " << pivot.str());
    random_rotation(pivot);
    perturb_rotate_->rotate_selection(pivot, rot_mat_, system);
  }

  virtual ~TrialPivot() {}

 private:
  int max_length_;
};

inline std::shared_ptr<TrialPivot> MakeTrialPivot(
    const argtype &args = argtype()) {
  return std::make_shared<TrialPivot>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_PIVOT_H_

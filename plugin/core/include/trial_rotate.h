
#ifndef FEASST_CORE_TRIAL_ROTATE_H_
#define FEASST_CORE_TRIAL_ROTATE_H_

#include "core/include/trial_move.h"
#include "core/include/perturb_rotate.h"
#include "core/include/random.h"
#include "core/include/utils_io.h"

namespace feasst {

class TrialRotate : public TrialMove {
 public:
  TrialRotate() {
    perturb_rotate_ = std::make_shared<PerturbRotate>();
    set_perturb(perturb_rotate_);
    set_recenter();
    set_tunable_param_max(180.);
    set_tunable_param_min(0.);
  }

  // randomize rotation matrix
  void random_rotation(
      ///pivot merely sets the dimensionality
      const Position& pivot) {
    Position axis = pivot;
    random_.unit_sphere_surface(&axis);
    const double angle = random_.uniform_real(-max_move(), max_move());
    rot_mat_.axis_angle(axis, angle);
  }

  void move(System * system) override {
    const Position& pivot = perturb_rotate_->selection().particle_positions()[0];
    random_rotation(pivot);
    perturb_rotate_->rotate_selection(pivot, rot_mat_, system);
  }

  void set_recenter(const double recenter = 0) { perturb_rotate_->set_recenter(recenter); }

  virtual ~TrialRotate() {}

 protected:
  // HWH they say not to use protected member variables but this makes
  // it easier to implement CrankShaft
  std::shared_ptr<PerturbRotate> perturb_rotate_;
  Random random_;
  RotationMatrix rot_mat_;
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_ROTATE_H_

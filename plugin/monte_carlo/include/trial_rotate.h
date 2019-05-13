
#ifndef FEASST_MONTE_CARLO_TRIAL_ROTATE_H_
#define FEASST_MONTE_CARLO_TRIAL_ROTATE_H_

#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "math/include/random.h"
#include "utils/include/utils_io.h"

namespace feasst {

class TrialRotate : public TrialMove {
 public:
  TrialRotate(
    /**
      max_move : for a given trial, the maximum possible size of the move.
     */
    const argtype& args = argtype()) : TrialMove(args) {
    perturb_rotate_ = std::make_shared<PerturbRotate>();
    // parse max move
    parse_tunable_(args, &args_, perturb_rotate_);
    perturb_rotate_->set_tune_min_and_max(0., 180.);
    set_perturb(perturb_rotate_);
    set_recenter();
  }

  void set_max_move_bounds(const Domain& domain) override {}

  // randomize rotation matrix
  void random_rotation(
      ///pivot merely sets the dimensionality
      const Position& pivot) {
    Position axis = pivot;
    random_.unit_sphere_surface(&axis);
    const double max_angle = perturb_rotate_->tunable().value();
    ASSERT(max_angle > NEAR_ZERO, "max angle is too small");
    const double angle = random_.uniform_real(-max_angle, max_angle);
    rot_mat_.axis_angle(axis, angle);
  }

  void move(Criteria * criteria, System * system) override {
    DEBUG("rotating");
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

inline std::shared_ptr<TrialRotate> MakeTrialRotate(
    const argtype &args = argtype()) {
  return std::make_shared<TrialRotate>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_ROTATE_H_

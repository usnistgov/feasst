
#ifndef FEASST_CORE_TRIAL_TRANSLATE_H_
#define FEASST_CORE_TRIAL_TRANSLATE_H_

#include "core/include/trial.h"
#include "core/include/perturb_translate.h"
#include "core/include/random.h"
#include "core/include/utils_io.h"

namespace feasst {

/**
 */
class TrialMove : public Trial {
 public:
  TrialMove() {
    set_group_index();
    // set_max_move();
  }
  void set_max_move(const double max_move = 0.1) {
    set_tunable_param(max_move); }
  double max_move() const { return tunable_param(); }
  void set_max_move_bounds(const Domain& domain) {
    set_tunable_param_max(domain.min_side_length()/2.);
    set_tunable_param_min(0.);
  }
  void set_perturb(std::shared_ptr<Perturb> perturb) { perturb_ = perturb; }
  virtual void move(System * system) = 0;
  virtual void select(System * system) {
    perturb_->select_random_particle(group_index(), system->configuration());
  }
  void attempt(Criteria* criteria, System * system) {
    before_attempt(criteria, system, perturb_.get());
    select(system);
    if (perturb_->selection().is_empty()) {
      // no particles present
      accept_criteria_.force_rejection = 1;
    } else {
      const double pe_old = system->energy(perturb_->selection());
      DEBUG("pe_old " << pe_old);
      move(system);
      const double pe_new = system->energy(perturb_->selection());
      DEBUG("pe_new " << pe_new);
      const double delta_energy = pe_new - pe_old;
      accept_criteria_.ln_metropolis_prob = -criteria->beta()*delta_energy;
      accept_criteria_.energy_new = criteria->running_energy() + delta_energy;
      accept_criteria_.energy_new_select = pe_new;
      accept_criteria_.force_rejection = 0;
      accept_criteria_.system = system;
      DEBUG("delta_energy " << delta_energy);
    }
    accept_or_reject(accept_criteria_, perturb_.get(), criteria);
  }
  virtual ~TrialMove() {}

 private:
  std::shared_ptr<Perturb> perturb_;
  AcceptanceCriteria accept_criteria_;
};

class TrialTranslate : public TrialMove {
 public:
  TrialTranslate() {
    perturb_trans_ = std::make_shared<PerturbTranslate>();
    set_perturb(perturb_trans_);
  }

  void move(System * system) override {
    const Position trajectory = random_.position_in_cube(
      system->dimension(),
      max_move()
    );
    perturb_trans_->translate_selection(trajectory, system);
  }

  virtual ~TrialTranslate() {}

 private:
  std::shared_ptr<PerturbTranslate> perturb_trans_;
  Random random_;
};

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

class TrialCrankShaft : public TrialRotate {
 public:
  TrialCrankShaft() : TrialRotate() { set_recenter(1); }

  void select(System * system) override {
    perturb_rotate_->select_random_segment_in_particle(group_index(), system->configuration());
  }

  void move(System * system) override {
    const Position& pivot = perturb_rotate_->selection().site_positions()[0].front();
    {
      Position axis = perturb_rotate_->selection().site_positions()[0].back();
      axis.subtract(pivot);
      axis.normalize();
      const double angle = random_.uniform_real(-max_move(), max_move());
      rot_mat_.axis_angle(axis, angle);
    }
    perturb_rotate_->rotate_selection(pivot, rot_mat_, system);
  }

  virtual ~TrialCrankShaft() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_TRANSLATE_H_

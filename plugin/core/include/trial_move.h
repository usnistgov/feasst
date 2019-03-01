
#ifndef FEASST_CORE_TRIAL_MOVE_H_
#define FEASST_CORE_TRIAL_MOVE_H_

#include "core/include/trial.h"
#include "core/include/perturb.h"
#include "core/include/utils_io.h"

namespace feasst {

/**
 */
class TrialMove : public Trial {
 public:
  TrialMove(
    /**
      max_move : for a given trial, the maximum possible size of the move.
     */
    const argtype &args = argtype()) : Trial(args) {
    // default
    set_group_index();
    // set_max_move();

    // parse
    args_.init(args);
    if (args_.key("max_move").used()) {
      set_max_move(args_.dble());
    }
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
  void before_attempt(Criteria* criteria, System * system, Perturb * perturb) override {
    Trial::before_attempt(criteria, system, perturb);
    set_max_move_bounds(system->configuration().domain());
  }
  void attempt(Criteria* criteria, System * system) override {
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

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_MOVE_H_

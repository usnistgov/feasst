
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
    const argtype &args = argtype()) : Trial(args) {
    // default
    set_group_index();
  }
  virtual void set_max_move_bounds(const Domain& domain) {};
  void set_perturb(std::shared_ptr<Perturb> perturb) { perturb_ = perturb; }
  virtual void move(Criteria * criteria, System * system) { ERROR("not implemented"); }
  virtual void move_and_acceptance(Criteria * criteria, System * system, AcceptanceCriteria * accept_criteria) {
    const double pe_old = system->energy(perturb_->selection());
    DEBUG("pe_old " << pe_old);
    move(criteria, system);
    if (accept_criteria_.force_rejection != 1) {
      const double pe_new = system->energy(perturb_->selection());
      DEBUG("pe_new " << pe_new);
      const double delta_energy = pe_new - pe_old;
      accept_criteria_.ln_metropolis_prob += -criteria->beta()*delta_energy;
      accept_criteria_.energy_new = criteria->running_energy() + delta_energy;
      accept_criteria_.energy_new_select = pe_new;
      accept_criteria_.system = system;
      DEBUG("delta_energy " << delta_energy);
    }
  }
  virtual void select(System * system) {
    perturb_->select_random_particle(group_index(), system->configuration());
  }
  void before_attempt(Criteria* criteria, System * system, Perturb * perturb, AcceptanceCriteria * accept_criteria) override {
    Trial::before_attempt(criteria, system, perturb, &accept_criteria_);
    set_max_move_bounds(system->configuration().domain());
  }
  void attempt(Criteria* criteria, System * system) override {
    before_attempt(criteria, system, perturb_.get(), &accept_criteria_);
    select(system);
    if (perturb_->selection().is_empty()) {
      // no particles present
      accept_criteria_.force_rejection = 1;
    } else {
      move_and_acceptance(criteria, system, &accept_criteria_);
    }
    accept_or_reject(accept_criteria_, perturb_.get(), criteria);
  }

  void tune() override {
    perturb_->tune(acceptance());
    reset_stats();
  }

  std::string status_header() const override {
    std::stringstream ss;
    ss << Trial::status_header() << " max_move";
    return ss.str();
  }

  std::string status() const override {
    std::stringstream ss;
    ss << Trial::status() << " " << perturb_->tunable().value();
    return ss.str();
  }

  virtual ~TrialMove() {}

  const std::shared_ptr<Perturb> perturb() const { return perturb_; }

 protected:
  void parse_tunable_(
      const argtype& args,
      Arguments * arguments,
      std::shared_ptr<Perturb> perturb) {
    arguments->init(args);
    Tunable tunable;
    tunable.set_value(arguments->key("max_move").dble());
    // tunable.set_value(arguments->key("max_move").dflt("0.1").dble());
    //if (arguments->key("max_move").used()) {
    //  tunable.set_value(arguments->dble());
    //}
    perturb->set_tunable(tunable);
  }

 private:
  std::shared_ptr<Perturb> perturb_;
  AcceptanceCriteria accept_criteria_;
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_MOVE_H_

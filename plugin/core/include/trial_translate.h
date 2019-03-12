
#ifndef FEASST_CORE_TRIAL_TRANSLATE_H_
#define FEASST_CORE_TRIAL_TRANSLATE_H_

#include "core/include/trial_move.h"
#include "core/include/perturb_translate.h"
#include "core/include/random.h"
#include "core/include/utils_io.h"
#include "core/include/rosenbluth.h"

namespace feasst {

class TrialTranslate : public TrialMove {
 public:
  TrialTranslate(
    /**
      max_move : for a given trial, the maximum possible size of the move.
     */
    const argtype& args = argtype()) : TrialMove(args) {
    translate_ = std::make_shared<PerturbTranslate>();
    // parse max move
    parse_tunable_(args, &args_, translate_);
    set_perturb(translate_);
  }

  void set_max_move_bounds(const Domain& domain) override {
    double min_length = domain.min_side_length()/2.;
    translate_->set_tune_min_and_max(0., min_length);
  }

  /// HWH move this to Trial
  void move(Criteria * criteria, System * system) override {
    translate_->perturb(system);
  }

  virtual ~TrialTranslate() {}

 private:
  std::shared_ptr<PerturbTranslate> translate_;
  Random random_;
};

inline std::shared_ptr<TrialTranslate> MakeTrialTranslate(const argtype &args = argtype()) {
  return std::make_shared<TrialTranslate>(args);
}

class StagedTrial {
 public:
  void parse_ref_and_num_steps(
      const argtype& args,
      Arguments * arguments,
      std::shared_ptr<Stage> stage) {
    arguments->init(args);
    if (arguments->key("reference").used()) {
      stage->set_reference(arguments->integer());
    }
    stage->set_num_steps(arguments->key("num_steps").integer());
  }
};

class TrialStagedTranslate : public TrialMove, public StagedTrial {
 public:
  TrialStagedTranslate(
    /**
      max_move : for a given trial, the maximum possible size of the move.
      reference : index of the reference potential.
      num_steps : number of steps per stage.
     */
    const argtype &args = argtype()) : TrialMove(args) {
    translate_ = std::make_shared<PerturbTranslate>();
    set_perturb(translate_);
    args_.init(args);
    parse_tunable_(args, &args_, translate_);
    auto stage = std::make_shared<Stage>();
    parse_ref_and_num_steps(args, &args_, stage);
    stage->set(translate_);
    // stage->set(std::make_shared<PerturbTranslate>(*translate_));
    stages_.add(stage);
  }

  void set_max_move_bounds(const Domain& domain) override {
    double min_length = domain.min_side_length()/2.;
    translate_->set_tune_min_and_max(0., min_length);
  }

  void select(System * system) override {
    TrialMove::select(system);
    stages_.parse_select(translate_->selection());
  }

  void move_and_acceptance(Criteria * criteria, System * system, AcceptanceCriteria * accept_criteria_) override {
    DEBUG("energy of select " << system->energy(translate_->selection()));
    int reject = 0;
    double en_old;
    const double ln_met_old = stages_.compute_rosenbluth(
      1,  // old
      translate_->selection(),
      criteria,
      system,
      &en_old,
      &reject);
    double en_new;
    const double ln_met_new = stages_.compute_rosenbluth(
      0,  // new
      translate_->selection(),
      criteria,
      system,
      &en_new,
      &reject);
    const double delta_energy = en_new - en_old;
    accept_criteria_->force_rejection = reject;
    accept_criteria_->ln_metropolis_prob -= ln_met_old;
    accept_criteria_->ln_metropolis_prob += ln_met_new;
    accept_criteria_->energy_new = criteria->running_energy() + delta_energy;
    accept_criteria_->energy_new_select = en_new;
    accept_criteria_->system = system;
    DEBUG("new en " << accept_criteria_->energy_new);
    DEBUG("delta_energy " << delta_energy);
    DEBUG("ln_met " << accept_criteria_->ln_metropolis_prob);
  }

  virtual ~TrialStagedTranslate() {}

 private:
  StageFactory stages_;
  std::shared_ptr<PerturbTranslate> translate_;
  Random random_;
};

inline std::shared_ptr<TrialStagedTranslate> MakeTrialStagedTranslate(
    const argtype &args = argtype()) {
  return std::make_shared<TrialStagedTranslate>(args);
}

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_TRANSLATE_H_

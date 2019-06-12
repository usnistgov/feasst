
#ifndef FEASST_MONTE_CARLO_TRIAL_H_
#define FEASST_MONTE_CARLO_TRIAL_H_

#include <vector>
#include <numeric>
#include <string>
#include <memory>
#include "monte_carlo/include/perturb.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/rosenbluth.h"

namespace feasst {

/**
  A stage contains both a selection and perturbation.
  Random realizations of the same perturbation is performed in a number of
  steps.
  A reference potential may be used during these steps, following the dual-
  cut configurational bias approach: http://doi.org/10.1080/002689798167881
 */
class TrialStage {
 public:
  TrialStage(
    /**
      num_steps: number of rosenbluth steps (default: 1).

      reference_index: index of reference potential.
        Otherwise, if full potential is desired, set to -1 (default: -1).
     */
    const argtype& args = argtype()) {
    Arguments args_(args);
    args_.dont_check();
    rosenbluth_.resize(args_.key("num_steps").dflt("1").integer());
    reference_ = args_.key("reference_index").dflt("-1").integer();
    set_mayer();
  }

  /// Return the index of the reference potential.
  int reference() const { return reference_; }

  /// Return true if the trial is utilizing Mayer sampling.
  bool is_mayer() const { return is_mayer_; }

  /// Set the above.
  void set_mayer(const bool enabled = false) { is_mayer_ = enabled; }

  /// Return the Rosenbluth.
  const Rosenbluth& rosenbluth() const { return rosenbluth_; }

  /// Set the selection.
  void set(std::shared_ptr<TrialSelect> select) { select_ = select; }

  /// Return the above.
  const TrialSelect * trial_select() const { return select_.get(); }

  /// Perform the selection.
  void select(System * system, Acceptance * acceptance) {
    select_->select(system);
    select_->check(acceptance);
  }

  /// Set the perturbation.
  void set(std::shared_ptr<Perturb> perturb) { perturb_ = perturb; }

//  /// Return the above.
//  const Perturb * perturb() const { return perturb_.get(); }

  /// Initialization before any stage attempt.
  void precompute(System * system) {
    select_->precompute(system);
    perturb_->precompute(select_.get(), system);
  }

  /// Initializations before each stage attempt.
  void before_stage() {
    perturb_->before_perturb();
  }

  /// Attempt all steps in a stage.
  /// Consider reference potentials and compute Rosenbluth factors.
  virtual void attempt(
    System * system,
    Criteria * criteria,
    /// Set to 1 for "old" system and "0" for new.
    const int old) {
    ASSERT(perturb_, "perturb not set");
    for (int step = 0; step < rosenbluth_.num(); ++step) {
      // DEBUG(perturb_->class_name());
      bool is_position_held = false;
      if (step == 0 and old == 1) is_position_held = true;
      perturb_->perturb(system, select_.get(), is_position_held);
      rosenbluth_.store(step, select_->mobile(), system);
      if (reference_ == -1) {
        DEBUG("select " << select_->mobile().str());
        rosenbluth_.set_energy(step, system->energy(select_->mobile()));
      } else {
        rosenbluth_.set_energy(step,
          system->reference_energy(select_->mobile(), reference_));
      }
      perturb_->revert(system);
    }
    rosenbluth_.compute(criteria->beta());
    if (old != 1) {
      if (rosenbluth_.chosen_step() != -1) {
        DEBUG("updating positions " << rosenbluth_.chosen().str());
        DEBUG("pos0 " << rosenbluth_.chosen().site_positions()[0][0].str());
        // DEBUG("pos1 " << rosenbluth_.chosen().site_positions()[0][1].str());
        system->get_configuration()->update_positions(rosenbluth_.chosen());
      } else if (is_mayer()) {
        ASSERT(rosenbluth_.num() == 1, "assumes 1 step for mayer");
        system->get_configuration()->update_positions(rosenbluth_.stored(0));
      }
    }
  }

  /// Revert the attempt.
  void revert(System * system) {
    perturb_->revert(system);
  }

  /// Finalize the attempt.
  void finalize(System * system) {
    DEBUG("finalizing");
    perturb_->finalize(system);
  }

  void tune(const double acceptance) {
    perturb_->tune(acceptance);
  }

 private:
  Rosenbluth rosenbluth_;
  int reference_ = -1;
  std::shared_ptr<Perturb> perturb_;
  std::shared_ptr<TrialSelect> select_;
  bool is_mayer_;
};

/// Implement the perturbation and calculation of acceptance.
class TrialCompute {
 public:
  /// Perform the stages on the system and compute the acceptance.
  void compute_rosenbluth(
    /// Set to 1 for "old" system and "0" for new.
    const int old,
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages) {
    double ln_rosenbluth = 0.;
    double energy_change = 0.;
    bool reference_used = false;
    for (TrialStage* stage : *stages) {
      stage->attempt(system, criteria, old);
      if (stage->rosenbluth().chosen_step() == -1) {
        if (!stage->is_mayer()) {
          acceptance->set_reject(true);
          DEBUG("auto reject");
          return;
        }
      }
      ln_rosenbluth += log(stage->rosenbluth().total_rosenbluth());
      double energy;
      if (old == 1) {
        energy = stage->rosenbluth().energy(0);
        acceptance->add_to_energy_old(energy);
      } else {
        energy = stage->rosenbluth().chosen_energy();
        acceptance->add_to_energy_new(energy);
      }
      energy_change += energy;
      if (stage->reference() >= 0) {
        reference_used = true;
      }
    }
    if (reference_used) {
      ASSERT(static_cast<int>(stages->size()) == 1,
        "implement a select_all in trial including selection of all trials?");
      const double en_full = system->energy((*stages)[0]->trial_select()->mobile());
      acceptance->add_to_ln_metropolis_prob(-1.*criteria->beta()*
        (en_full - energy_change));
      DEBUG("energy ref: " << energy_change);
      acceptance->set_energy_ref(energy_change);
      if (old == 1) {
        acceptance->set_energy_old(en_full);
      } else {
        acceptance->set_energy_new(en_full);
      }
    }
    acceptance->add_to_ln_metropolis_prob(ln_rosenbluth);
  }

  virtual void perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages) = 0;
};

/*
  Trial optionals:
    - internal vs external potential (e.g., bonded cb, long-range, etc)

  Trial:
    - before trial attempts: precompute
    - during trial attempts:
    - select all
      - turn off interactions for selections which have not been perturbed yet
      - add can't select anything before perturbation
    - cycle CB stages - not, CB may require old/new building
      - select stage
      - record old (opt: ref)
      - cycle k steps
        - perturb(select_stage)
        - record new (opt: ref)
        - if k!=1, revert
        - if k==1, continue perturbations (even multiple stages before decision?)
      - update rosenbluth
      - turn on interaction of perturbed stage selection
      - if k!=1,
        - decide
        - replicate chosen perturb
    - decide all stages (opt: full - ref)
      - note, mayer sampling can use ref term computed during stages
      - note, current energy.. may be obtained in different ways
        - for most, current = previous current + delta
        - for mayer sampling, don't use delta energy.
    - keep/revert all stages

 */

/**
  A trial contains a number of TrialStages.
  The Acceptance is computed as the stages are enacted, and then sent to
  Criteria to decide if the trial is accepted or rejected.
 */
class Trial {
 public:
  Trial() {
    set_weight();
    set_mayer();
  }

  /// Return the unnormalized relative probability of selection of this trial
  /// with respect to all trials.
  double weight() const { return weight_; }

  /// Set the above weight.
  void set_weight(const double weight = 1.) { weight_ = weight; }

  /// Add a TrialStage in order.
  // HWH depreciate
  void add(std::shared_ptr<TrialStage> stage) {
    stages_.push_back(stage);
    stages_ptr_.push_back(stages_.back().get());
  }

  /// Add a stage which includes selection and perturbation with arguments.
  void add_stage(
    std::shared_ptr<TrialSelect> select,
    std::shared_ptr<Perturb> perturb,
    const argtype& args = argtype()) {
    auto stage = std::make_shared<TrialStage>(args);
    stage->set(select);
    stage->set(perturb);
    add(stage);
  }

  /// Return the stages.
  const std::vector<std::shared_ptr<TrialStage> > stages() const { return stages_; }

  /// Return a stage.
  const TrialStage * stage(const int index) const { return stages_[index].get(); }

  /// Number of stages.
  int num_stages() const { return static_cast<int>(stages_.size()); }

  /// Number of successful attempts.
  int64_t num_success() const { return num_success_; }

  /// Number of attempts.
  int64_t num_attempts() const { return num_attempts_; }

  void increment_num_attempts() { ++num_attempts_; }

  /// Return the ratio of the number of successful attempts and total attempts.
  double acceptance() const {
    return static_cast<double>(num_success_)/static_cast<double>(num_attempts_);
  }

  /// Reset trial statistics.
  virtual void reset_stats() {
    num_attempts_ = 0;
    num_success_ = 0;
  }

  /// Return the header description for the status of the trial (e.g., acceptance, etc).
  virtual std::string status_header() const {
    std::stringstream ss;
    ss << "acceptance";
    return ss.str();
  }

  /// Return the status of the trial (e.g., acceptance, etc).
  virtual std::string status() const {
    std::stringstream ss;
    ss << acceptance();
    return ss.str();
  }

  virtual void tune() {
    for (auto stage : stages_) stage->tune(acceptance());
    reset_stats();
  }

  /// Set a Mayer-sampling trial.
  void set_mayer(const bool enabled = false) {
    for (auto stage : stages_) stage->set_mayer(enabled); }

  /// Precompute quantities before simulation for optimization.
  virtual void precompute(Criteria * criteria, System * system) {
    for (std::shared_ptr<TrialStage> stage : stages_) {
      stage->precompute(system);
    }
  }

  /// Set the computation of the trial and acceptance.
  void set(std::shared_ptr<TrialCompute> compute) { compute_ = compute; }

  /// Attempt a trial.
  virtual void attempt(Criteria * criteria, System * system) {
    ++num_attempts_;
    acceptance_.reset();
    criteria->before_attempt(system);
    for (std::shared_ptr<TrialStage> stage : stages_) {
      stage->before_stage();
      stage->select(system, &acceptance_);
    }
    if (!acceptance_.reject()) {
      compute_->perturb_and_acceptance(
        criteria, system, &acceptance_, &stages_ptr_);
    }
    if (criteria->is_accepted(acceptance_, system)) {
      DEBUG("accepted");
      ++num_success_;
      for (int index = num_stages() - 1; index >= 0; --index) {
        stages_[index]->finalize(system);
      }
    } else {
      DEBUG("rejected");
      for (int index = num_stages() - 1; index >= 0; --index) {
        stages_[index]->revert(system);
      }
    }
  }

  std::string class_name() const { return class_name_; }

 protected:
  std::string class_name_ = "Trial";

 private:
  std::vector<std::shared_ptr<TrialStage> > stages_;
  std::vector<TrialStage*> stages_ptr_;
  std::shared_ptr<TrialCompute> compute_;
  Acceptance acceptance_;
  double weight_ = 1.;
  int64_t num_attempts_ = 0, num_success_ = 0;
};

class TrialComputeMove : public TrialCompute {
 public:
  void perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages) override {
    compute_rosenbluth(1, criteria, system, acceptance, stages);
    compute_rosenbluth(0, criteria, system, acceptance, stages);
    acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new() - acceptance->energy_old());
  }
};

/// Attempt to rigidly move a selection in a Trial in one stage.
class TrialMove : public Trial {
 public:
  TrialMove(
    std::shared_ptr<TrialSelect> select,
    std::shared_ptr<PerturbMove> perturb,
    const argtype& args = argtype()) {
    add_stage(select, perturb, args);
    set(std::make_shared<TrialComputeMove>());
  }
};

/// Attempt a rigid translation of a random particle.
class TrialTranslate : public TrialMove {
 public:
  TrialTranslate(
    /// These arguments are sent to both PerturbTranslate and TrialStage.
    const argtype& args = argtype())
    : TrialMove(
      std::make_shared<TrialSelectParticle>(),
      std::make_shared<PerturbTranslate>(args),
      args
    ) {};
};

inline std::shared_ptr<TrialTranslate> MakeTrialTranslate(
    const argtype &args = argtype()) {
  return std::make_shared<TrialTranslate>(args);
}

/// Attempt a rigid rotation of a random particle.
class TrialRotate : public TrialMove {
 public:
  TrialRotate(
    /// These arguments are sent to both PerturbRotate and TrialStage.
    const argtype& args = argtype())
    : TrialMove(
      std::make_shared<TrialSelectParticle>(),
      std::make_shared<PerturbRotate>(args),
      args
    ) {};
};

inline std::shared_ptr<TrialRotate> MakeTrialRotate(
    const argtype &args = argtype()) {
  return std::make_shared<TrialRotate>(args);
}

class TrialComputeAdd : public TrialCompute {
 public:
  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages) override {
    compute_rosenbluth(0, criteria, system, acceptance, stages);
    const TrialSelect * select = (*stages)[0]->trial_select();
    system->get_configuration()->revive(select->mobile());
    acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
    { // Metropolis
      const Configuration& config = system->configuration();
      const double volume = config.domain().volume();
      const int particle_type = select->particle_type();
      acceptance->add_to_ln_metropolis_prob(
        log(volume*select->probability())
        + criteria->beta_mu(particle_type)
      );
    }
  }
};

/// Attempt to add a particle.
class TrialAdd : public Trial {
 public:
  TrialAdd(const argtype& args = argtype()) {
    auto perturb = std::make_shared<PerturbAdd>(args);
    add_stage(
      std::make_shared<TrialSelectDoNothing>(),
      perturb,
      args
    );
    set(std::make_shared<TrialComputeAdd>());
  }
};

inline std::shared_ptr<TrialAdd> MakeTrialAdd(
    const argtype &args = argtype()) {
  return std::make_shared<TrialAdd>(args);
}

class TrialComputeRemove : public TrialCompute {
 public:
  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages) override {
    compute_rosenbluth(1, criteria, system, acceptance, stages);
    acceptance->set_energy_new(criteria->current_energy() - acceptance->energy_old());
    acceptance->set_macrostate_shift(-1);
    { // Metropolis
      const Configuration& config = system->configuration();
      const double volume = config.domain().volume();
      const TrialSelect * select = (*stages)[0]->trial_select();
      const int particle_type = select->particle_type();
      DEBUG("volume " << volume << " selprob " << select->probability() << " betamu " << criteria->beta_mu(particle_type));
      acceptance->add_to_ln_metropolis_prob(
        - log(volume*select->probability())
        - criteria->beta_mu(particle_type)
      );
      DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    }
  }
};

/// Attempt to remove a particle.
class TrialRemove : public Trial {
 public:
  TrialRemove(const argtype& args = argtype()) {
    argtype args2(args);
    args2.insert({"load_coordinates", "false"});
    add_stage(
      std::make_shared<TrialSelectParticleOfType>(args2),
      std::make_shared<PerturbRemove>(),
      args
    );
    set(std::make_shared<TrialComputeRemove>());
  }
};

inline std::shared_ptr<TrialRemove> MakeTrialRemove(
    const argtype &args = argtype()) {
  return std::make_shared<TrialRemove>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_H_

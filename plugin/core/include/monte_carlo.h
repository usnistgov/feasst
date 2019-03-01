
#ifndef FEASST_CORE_MONTE_CARLO_H_
#define FEASST_CORE_MONTE_CARLO_H_

#include <vector>
#include <memory>
#include "core/include/trial_factory.h"
#include "core/include/trial_transfer.h"
#include "core/include/analyze.h"
#include "core/include/modify.h"

namespace feasst {

/**
  Enforce order of Monte Carlo initialization:
  1 config
    - particle types
    - domain
    - model parameters (based on types)
    - groups
  2 potentials
    - models
    - model precomputes(type depend)
    - visitors
      - cells
    - visitor precomputes
  3 system finalized
  4 criteria
    - set running energy
  5 trials and steppers in mixed order
    - trial precomputes?
    - nmolseek
    - analyzers and modifers
      - ana and mod precomputes
 */
class MonteCarlo {
 public:
  /// The first action with a Monte Carlo object is to set the System.
  /// This must be done before setting Criteria.
  void set(const System& system) {
    system_set_ = true;
    system_ = system;
  }

  /// Once the System is set, it may be accessed on a read-only basis.
  const System& system() const { return system_; }

  // HWH depreciate: only in rare cases should the system be modified directly.
  System * get_system() { return &system_; }

  /// The second action is to set the Criteria.
  /// System must be set first.
  void set(std::shared_ptr<Criteria> criteria) {
    ASSERT(system_set_, "set System before Criteria.");
    criteria->set_running_energy(system_.energy());
    criteria_ = criteria;
    criteria_set_ = true;
  }

  /// Once Criteria is set, it may be accessed on a read-only basis.
  const std::shared_ptr<Criteria> criteria() { return criteria_; }

  /// The remaining actions can be done in almost any order.
  /// Typically, one begins by adding trials.
  void add(std::shared_ptr<Trial> trial) {
    ASSERT(criteria_set_, "set Criteria before Trials.");
    trial_factory_.add(trial);
    // If later, perhaps after some initialization, more trials are added,
    // then Analyze and Modify classes may need to be re-initialized.
    // analyze_factory_.initialize(criteria_, system_, trial_factory_);
    // modify_factory_.initialize(criteria_, &system_, &trial_factory_);
  }

  /// Access the trials on a read-only basis.
  const TrialFactory& trials() const { return trial_factory_; }

  /// HWH consider depreciating this. While convenient, it confuses many users.
  /// Before the production simulation, one may wish to quickly seek an initial
  /// configuration with a given number of particles, without necessarily
  /// satisfying detailed balance.
  void seek_num_particles(const int num) {
    TrialTransfer add;
    add.set_add_probability(1);
    while (system_.configuration().num_particles() < num) {
      attempt();
      add.attempt(criteria_.get(), &system_);
    }
    trial_factory_.reset_stats();
  }

  /// An Analyzer performs some task after a given number of steps, but is
  /// read-only on System, Criteria and Trials.
  void add(const std::shared_ptr<Analyze> analyze) {
    analyze->initialize(criteria_, system_, trial_factory_);
    analyze_factory_.add(analyze);
  }

  /// A Modifier performs some task after a given number of steps, but may
  /// change the System, Criteria and Trials.
  void add(const std::shared_ptr<Modify> modify) {
    modify->initialize(criteria_, &system_, &trial_factory_);
    modify_factory_.add(modify);
  }

  /// Attempt a number of Monte Carlo trials, with subsequent Analyzers and
  /// Modifiers.
  void attempt(const int num_trials = 1) {
    for (int trial = 0; trial < num_trials; ++trial) {
      trial_factory_.attempt(criteria_.get(), &system_);
      analyze_factory_.trial(criteria_, system_, trial_factory_);
      modify_factory_.trial(criteria_, &system_, &trial_factory_);
    }
  }

 private:
  std::shared_ptr<Criteria> criteria_;
  TrialFactory trial_factory_;
  System system_;
  AnalyzeFactory analyze_factory_;
  ModifyFactory modify_factory_;

  bool system_set_ = false;
  bool criteria_set_ = false;
};

}  // namespace feasst

#endif  // FEASST_CORE_MONTE_CARLO_H_

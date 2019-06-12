
#ifndef FEASST_MONTE_CARLO_MONTE_CARLO_H_
#define FEASST_MONTE_CARLO_MONTE_CARLO_H_

#include <sstream>
#include <iostream>
#include <vector>
#include <memory>
#include "utils/include/checkpoint.h"
//#include "utils/include/timer.h"
#include "monte_carlo/include/trial_factory.h"
//#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/analyze.h"
#include "monte_carlo/include/analyze_factory.h"
#include "monte_carlo/include/modify.h"
#include "monte_carlo/include/modify_factory.h"

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
  MonteCarlo() {}
//    timer_other_ = timer_.add("other");
//    timer_trial_ = timer_.add("trial");
//    timer_analyze_ = timer_.add("analyze");
//    timer_modify_ = timer_.add("modify");
//    timer_checkpoint_ = timer_.add("checkpoint");
//  }

  /// The first action with a Monte Carlo object is to set the Configuration.
  void add(const Configuration& config) {
    system_.add(config);
    config_set_ = true;
    if (potential_set_) system_set_ = true;
    ASSERT(!criteria_set_, "add config before criteria");
  }

  /// The configuration may be accessed read-only.
  const Configuration& configuration(const int index = 0) const {
    return system_.configuration(index); }

  /// The second action is to add Potentials.
  void add(const Potential& potential) {
    ASSERT(!criteria_set_, "add potential before criteria");
    system_.add(potential);
    potential_set_ = true;
    if (config_set_) system_set_ = true;
  }

  /// Add potential to reference.
  void add_to_reference(const Potential& potential) {
    system_.add_to_reference(potential);
  }

  /// Alternatively, the first and second actions may be combined by setting
  /// the system directly.
  /// This must be done before setting Criteria.
  void set(const System& system) {
    system_set_ = true;
    system_ = system;
    system_.precompute();
    ASSERT(!criteria_set_, "add system before criteria");
  }

  /// Once the System is set, it may be accessed on a read-only basis.
  const System& system() const { return system_; }

  // HWH depreciate: only in rare cases should the system be modified directly.
  System * get_system() { return &system_; }
  Criteria * get_criteria() { return criteria_.get(); }

  /// The third action is to set the Criteria.
  /// Configuration and Potentials (or System) must be set first.
  void add(std::shared_ptr<Criteria> criteria) {
    ASSERT(system_set_, "set System before Criteria.");
    criteria->set_current_energy(system_.energy());
    criteria_ = criteria;
    criteria_set_ = true;
  }

  // HWH depreciate
  void set(std::shared_ptr<Criteria> criteria) { add(criteria); }

  /// Once Criteria is set, it may be accessed on a read-only basis.
  const std::shared_ptr<Criteria> criteria() { return criteria_; }

  /// The remaining actions can be done in almost any order.
  /// Typically, one begins by adding trials.
  void add(std::shared_ptr<Trial> trial) {
    ASSERT(criteria_set_, "set Criteria before Trials.");
    trial->precompute(criteria_.get(), &system_);
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
    ASSERT(system_.configuration().num_particles() <= num,
      "assumes you only want to add particles, not delete them");
    TrialAdd add;
    while (system_.configuration().num_particles() < num) {
      attempt();
      add.attempt(criteria_.get(), &system_);
    }
    trial_factory_.reset_stats();
  }

  /// An Analyzer performs some task after a given number of steps, but is
  /// read-only on System, Criteria and Trials.
  /// At this stage, multistate non-factory classes are converted into
  /// factories for each state in criteria.
  void add(std::shared_ptr<Analyze> analyze);

  std::vector<std::shared_ptr<Analyze> > analyzers() const {
    return analyze_factory_.analyzers(); }
  std::shared_ptr<Analyze> analyze(const int index) const {
    return analyze_factory_.analyzers()[index]; }

  /// A Modifier performs some task after a given number of steps, but may
  /// change the System, Criteria and Trials.
  void add(const std::shared_ptr<Modify> modify) {
    ASSERT(criteria_set_, "set Criteria before Modify");
    modify->initialize(criteria_, &system_, &trial_factory_);
    modify_factory_.add(modify);
  }

  /// Add a checkpoint.
  void add(const std::shared_ptr<Checkpoint> checkpoint) {
    checkpoint_ = checkpoint; }

  /// Attempt a number of Monte Carlo trials, with subsequent Analyzers and
  /// Modifiers.
  void attempt(const int num_trials = 1) {
    ASSERT(system_set_, "system must be set before attempting trials.");
    ASSERT(criteria_set_, "criteria must be set before attempting trials.");
    for (int trial = 0; trial < num_trials; ++trial) {
      DEBUG("mc trial: " << trial);
      //timer_.start(timer_trial_);
      trial_factory_.attempt(criteria_.get(), &system_);
      //timer_.start(timer_analyze_);
      analyze_factory_.trial(criteria_, system_, trial_factory_);
      //timer_.start(timer_modify_);
      modify_factory_.trial(criteria_, &system_, &trial_factory_);
      if (checkpoint_) {
        //timer_.start(timer_checkpoint_);
        checkpoint_->check(*this);
      }
    }
    //timer_.start(timer_other_);
  }

  /// Attempt Monte Carlo trials until Criteria returns completion.
  void run_until_complete() {
    while (!criteria_->is_complete()) {
      attempt(1);
    }
  }

//  const Timer& timer() const { return timer_; }
//  std::string timer_str() const {
//    std::stringstream ss;
//    const double trial_missing = timer_.missing_percent("trial", trial_factory_.timer());
//    ss << timer_.str()
//       << "*** TrialFactory Profile ***" << std::endl
//       << trial_factory_.timer().str()
//       << "missing CPU hours percentage: " << trial_missing;
//    return ss.str();
//  }

  void serialize(std::ostream& ostr) const {
    feasst_serialize_version(1, ostr);
    feasst_serialize_fstobj(system_, ostr);
    feasst_serialize_fstdr(criteria_, ostr);
    // feasst_serialize_fstobj(trial_factory_, ostr);
    feasst_serialize_fstobj(analyze_factory_, ostr);
    feasst_serialize_fstobj(modify_factory_, ostr);
    feasst_serialize(checkpoint_, ostr);
  }

  MonteCarlo(std::istream& istr) {
    feasst_deserialize_version(istr);
    feasst_deserialize_fstobj(&system_, istr);
    // feasst_deserialize_fstdr(criteria_, istr);
    { // HWH for unknown reasons the above template function does not work
      int existing;
      istr >> existing;
      if (existing != 0) {
        criteria_ = criteria_->deserialize(istr);
      }
    }
    // feasst_deserialize_fstobj(&trial_factory_, istr);
    feasst_deserialize_fstobj(&analyze_factory_, istr);
    feasst_deserialize_fstobj(&modify_factory_, istr);
    feasst_deserialize(checkpoint_, istr);
  }

 private:
  System system_;
  std::shared_ptr<Criteria> criteria_;
  TrialFactory trial_factory_;
  AnalyzeFactory analyze_factory_;
  ModifyFactory modify_factory_;
  std::shared_ptr<Checkpoint> checkpoint_;

//  Timer timer_;
//  int timer_other_, timer_trial_, timer_analyze_, timer_modify_;
//  int timer_checkpoint_;

  bool config_set_ = false;
  bool potential_set_ = false;
  bool system_set_ = false;
  bool criteria_set_ = false;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_MONTE_CARLO_H_

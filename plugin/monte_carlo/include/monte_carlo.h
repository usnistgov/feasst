
#ifndef FEASST_MONTE_CARLO_MONTE_CARLO_H_
#define FEASST_MONTE_CARLO_MONTE_CARLO_H_

#include <sstream>
#include <iostream>
#include <vector>
#include <memory>
#include "utils/include/checkpoint.h"
#include "math/include/random_mt19937.h"
//#include "utils/include/timer.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_remove.h"
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
  MonteCarlo() {
    random_ = std::make_shared<RandomMT19937>();
//    timer_other_ = timer_.add("other");
//    timer_trial_ = timer_.add("trial");
//    timer_analyze_ = timer_.add("analyze");
//    timer_modify_ = timer_.add("modify");
//    timer_checkpoint_ = timer_.add("checkpoint");
  }

  /// Set the random number generator.
  void set(std::shared_ptr<Random> random) { random_ = random; }

  /// Return the random number generator.
  const Random * random() const { return random_.get(); }

  /// Offset random number generator.
  void offset_random() { random_->uniform(); }

  /// The first action with a Monte Carlo object is to set the Configuration.
  void add(const Configuration& config);

  /// The configuration may be accessed read-only.
  const Configuration& configuration(const int index = 0) const {
    return system_.configuration(index); }

  /// The second action is to add Potentials.
  void add(const Potential& potential);

  /// Add potential to optimized.
  void add_to_optimized(const Potential& potential) {
    system_.add_to_optimized(potential); }

  /// Add potential to reference.
  void add_to_reference(const Potential& potential) {
    system_.add_to_reference(potential); }

  /// Alternatively, the first and second actions may be combined by setting
  /// the system directly.
  /// This must be done before setting Criteria.
  void set(const System& system);

  /// Once the System is set, it may be accessed on a read-only basis.
  const System& system() const { return system_; }

  // HWH depreciate: only in rare cases should the system be modified directly.
  System * get_system() { return &system_; }
  Criteria * get_criteria() { return criteria_.get(); }

  /// The third action is to set the Criteria.
  /// Configuration and Potentials (or System) must be set first.
  void add(std::shared_ptr<Criteria> criteria);

  // HWH depreciate
  void set(std::shared_ptr<Criteria> criteria) { add(criteria); }

  /// Once Criteria is set, it may be accessed on a read-only basis.
  const Criteria * criteria() const { return criteria_.get(); }

  /// The remaining actions can be done in almost any order.
  /// Typically, one begins by adding trials.
  void add(std::shared_ptr<Trial> trial);

  /// Access the trials on a read-only basis.
  const TrialFactory& trials() const { return trial_factory_; }

  /// Access the trials on a read-only basis.
  const Trial * trial(const int index) const {
    return trial_factory_.trial(index); }

  /// HWH consider depreciating this. While convenient, it confuses many users.
  /// Before the production simulation, one may wish to quickly seek an initial
  /// configuration with a given number of particles, without necessarily
  /// satisfying detailed balance.
  /// Thus, all stats are reset.
  void seek_num_particles(const int num);

  /// An Analyzer performs some task after a given number of steps, but is
  /// read-only on System, Criteria and Trials.
  /// At this stage, multistate non-factory classes are converted into
  /// factories for each state in criteria.
  void add(std::shared_ptr<Analyze> analyze);

  /// Return all analyzers.
  const std::vector<std::shared_ptr<Analyze> >& analyzers() const {
    return analyze_factory_.analyzers(); }

  /// Return an Analayze by index.
  const Analyze * analyze(const int index) const {
    return analyze_factory_.analyzers()[index].get(); }

  /// Return the number of analyzers.
  int num_analyzers() const {
    return static_cast<int>(analyze_factory_.analyzers().size()); }

  /// A Modifier performs some task after a given number of steps, but may
  /// change the System, Criteria and Trials.
  void add(const std::shared_ptr<Modify> modify);

  /// Add a checkpoint.
  void add(const std::shared_ptr<Checkpoint> checkpoint) {
    checkpoint_ = checkpoint; }

  void before_attempts_() {
    ASSERT(system_set_, "system must be set before attempting trials.");
    ASSERT(criteria_set_, "criteria must be set before attempting trials.");
  }

  /// Attempt one trial, with subsequent analysers and modifiers.
  void attempt() {
//    before_attempts_();
    attempt_(1, &trial_factory_, random_.get());
    //timer_.start(timer_trial_);
//    trial_factory_.attempt(criteria, system, random);
//    after_trial_();
  }

  /// Attempt a number of Monte Carlo trials.
  void attempt(const int num_trials) { attempt_(num_trials, &trial_factory_, random_.get()); }

  /// Reset trial statistics
  virtual void reset_trial_stats() { trial_factory_.reset_stats(); }

  /// Attempt trial index without analyzers, modifiers or checkpoints.
  bool attempt_trial(const int index) {
    return trial_factory_.attempt(criteria_.get(), &system_, index, random_.get());
  }

  /// Revert changes from previous trial.
  void revert(const int trial_index, const bool accepted) {
    trial_factory_.revert(trial_index, accepted, &system_);
    DEBUG("reverting " << criteria_->current_energy());
    criteria_->revert(accepted);
  }

  /// Finalize changes from previous trial.
  void finalize(const int trial_index) {
    trial_factory_.finalize(trial_index, &system_);
  }

  // HWH hackish interface for pipeline
  void delay_finalize() {
    trial_factory_.delay_finalize();
  }

  /// Attempt Monte Carlo trials until Criteria returns completion.
  void run_until_complete() {
    while (!criteria_->is_complete()) {
      attempt(1);
    }
  }

  /// Mimic a rejection by a trial.
  void mimic_trial_rejection(const int trial_index, const double ln_prob) {
    trial_factory_.mimic_trial_rejection(trial_index);
    // HWH add this for FH
    // criteria_->mimic_trial_rejection(ln_prob);
  }

  /// Load random numbers and energy calculations into cache.
  void load_cache(const bool load) {
    random_->set_cache_to_load(load);
    system_.load_cache(load);
  }

  /// Unload random numbers and energy calculations from cache.
  void unload_cache(const MonteCarlo& mc) {
  //Random& random, const System& system) {
    random_->set_cache_to_unload((*mc.random_));
    system_.unload_cache(mc.system());
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

  void after_trial_() {
    analyze_factory_.trial(criteria_.get(), system_, trial_factory_);
    modify_factory_.trial(criteria_.get(), &system_, &trial_factory_);
    if (checkpoint_) {
      checkpoint_->check(*this);
    }
  }

  virtual void serialize(std::ostream& ostr) const {
    feasst_serialize_version(529, ostr);
    feasst_serialize_fstobj(system_, ostr);
    feasst_serialize_fstdr(criteria_, ostr);
    feasst_serialize_fstobj(trial_factory_, ostr);
    feasst_serialize_fstobj(analyze_factory_, ostr);
    feasst_serialize_fstobj(modify_factory_, ostr);
    feasst_serialize(checkpoint_, ostr);
    feasst_serialize_fstdr(random_, ostr);
    feasst_serialize(config_set_, ostr);
    feasst_serialize(potential_set_, ostr);
    feasst_serialize(system_set_, ostr);
    feasst_serialize(criteria_set_, ostr);
  }

  std::string serialize() {
    std::stringstream ss;
    serialize(ss);
    return ss.str();
  }

  MonteCarlo(std::istream& istr) {
    //INFO(istr.rdbuf());
    //int tmp;
    //istr >> tmp;
    //INFO("tmp " << tmp);
    const int version = feasst_deserialize_version(istr);
    ASSERT(version == 529, "version: " << version);
    feasst_deserialize_fstobj(&system_, istr);
    // feasst_deserialize_fstdr(criteria_, istr);
    { // HWH for unknown reasons the above template function does not work
      int existing;
      istr >> existing;
      if (existing != 0) {
        criteria_ = criteria_->deserialize(istr);
      }
    }
    feasst_deserialize_fstobj(&trial_factory_, istr);
    feasst_deserialize_fstobj(&analyze_factory_, istr);
    feasst_deserialize_fstobj(&modify_factory_, istr);
    feasst_deserialize(checkpoint_, istr);
    // HWH for unknown reasons, this function template does not work.
    //feasst_deserialize_fstdr(random_, istr);
    { int existing;
      istr >> existing;
      if (existing != 0) {
        random_ = random_->deserialize(istr);
      }
    }
    feasst_deserialize(&config_set_, istr);
    feasst_deserialize(&potential_set_, istr);
    feasst_deserialize(&system_set_, istr);
    feasst_deserialize(&criteria_set_, istr);
  }

  MonteCarlo deserialize(const std::string str) {
    std::stringstream ss(str);
    return MonteCarlo(ss);
  }

  virtual ~MonteCarlo() {}

 protected:
  virtual void attempt_(const int num_trials, TrialFactory * trial_factory, Random * random) {
    before_attempts_();
    for (int trial = 0; trial < num_trials; ++trial) {
      DEBUG("mc trial: " << trial);
      trial_factory->attempt(criteria_.get(), &system_, random_.get());
      after_trial_();
    }
  }

 private:
  System system_;
  std::shared_ptr<Criteria> criteria_;
  TrialFactory trial_factory_;
  AnalyzeFactory analyze_factory_;
  ModifyFactory modify_factory_;
  std::shared_ptr<Checkpoint> checkpoint_;
  std::shared_ptr<Random> random_;

//  Timer timer_;
//  int timer_other_, timer_trial_, timer_analyze_, timer_modify_;
//  int timer_checkpoint_;

  bool config_set_ = false;
  bool potential_set_ = false;
  bool system_set_ = false;
  bool criteria_set_ = false;
};

/// Initialize an add and remove trial simultaneously with the same arguments
inline void add_trial_transfer(MonteCarlo * mc, const argtype& args = argtype()) {
  mc->add(MakeTrialAdd(args));
  mc->add(MakeTrialRemove(args));
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_MONTE_CARLO_H_

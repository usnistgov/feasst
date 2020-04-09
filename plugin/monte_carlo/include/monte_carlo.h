
#ifndef FEASST_MONTE_CARLO_MONTE_CARLO_H_
#define FEASST_MONTE_CARLO_MONTE_CARLO_H_

#include <sstream>
#include <vector>
#include <memory>
//#include "utils/include/timer.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/analyze.h"
#include "monte_carlo/include/analyze_factory.h"
#include "monte_carlo/include/modify.h"
#include "monte_carlo/include/modify_factory.h"

namespace feasst {

class Checkpoint;
class Random;
class RandomMT19937;

// HWH document this class better
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
  /// Construct a MonteCarlo object with RandomMT19937.
  MonteCarlo();

  /// Set the random number generator.
  void set(std::shared_ptr<Random> random) { random_ = random; }

  /// Return the random number generator.
  const Random * random() const { return random_.get(); }

  /// Seed random number generator.
  void seed_random(const int seed);

  /// The first action with a Monte Carlo object is to set the Configuration.
  void add(const Configuration& config);

  /// The configuration may be accessed read-only.
  const Configuration& configuration(const int index = 0) const {
    return system_.configuration(index); }

  /// The second action is to add Potentials.
  void add(const Potential& potential);

  /// Set an existing potential.
  void set(const int index, const Potential& potential);

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
  void add(const std::shared_ptr<Checkpoint> checkpoint);

  void before_attempts_();

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
  void revert(const int trial_index,
      const bool accepted,
      const double ln_prob);

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
    run_until_complete_(&trial_factory_, random_.get());
  }

  /// Mimic a rejection by a trial.
  void imitate_trial_rejection(const int trial_index,
      const double ln_prob,
      const int state_old,
      const int state_new) {
    trial_factory_.imitate_trial_rejection(trial_index);
    criteria_->imitate_trial_rejection(ln_prob, state_old, state_new);
  }

  /// Load random numbers and energy calculations into cache.
  void load_cache(const bool load);

  /// Unload random numbers and energy calculations from cache.
  void unload_cache(const MonteCarlo& mc);

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
    after_trial_analyze_();
    after_trial_modify_();
  }

  void after_trial_analyze_() {
    analyze_factory_.trial(criteria_.get(), system_, trial_factory_);
  }

  void after_trial_modify_();

  const PhysicalConstants * physical_constants() const {
    return configuration().model_params().physical_constants(); }

  virtual void serialize(std::ostream& ostr) const;

  std::string serialize() {
    std::stringstream ss;
    serialize(ss);
    return ss.str();
  }

  explicit MonteCarlo(std::istream& istr);

  MonteCarlo deserialize(const std::string str) {
    std::stringstream ss(str);
    return MonteCarlo(ss);
  }

  virtual ~MonteCarlo() {}

 protected:
  virtual void attempt_(int num_trials,
    TrialFactory * trial_factory,
    Random * random);

  virtual void run_until_complete_(TrialFactory * trial_factory,
                                   Random * random) {
    while (!criteria_->is_complete()) {
      attempt_(1, trial_factory, random);
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

/// Initialize an add and remove trial simultaneously with the same arguments.
void add_trial_transfer(MonteCarlo * mc, const argtype& args = argtype());

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_MONTE_CARLO_H_


#ifndef FEASST_MONTE_CARLO_MONTE_CARLO_H_
#define FEASST_MONTE_CARLO_MONTE_CARLO_H_

#include <sstream>
#include <vector>
#include <memory>
//#include "utils/include/timer.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/analyze.h"
#include "monte_carlo/include/analyze_factory.h"
#include "monte_carlo/include/modify.h"
#include "monte_carlo/include/modify_factory.h"

namespace feasst {

class Checkpoint;
class Random;

// HWH document this class better
// HWH consider a MonteCarlo interface where everything is input through
// HWH constructors, so that order of initialization is thoroughly enforced
// HWH requires ability to initialize simulation (e.g., seek_num)
// HWH or control of when Analyze/Modify begin (e.g., init/equil phases)
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
  //const Random& random() const { return const_cast<Random&>(*random_); }

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
  Random * get_random() { return random_.get(); }
  TrialFactory * get_trial_factory() { return &trial_factory_; }

  /// The third action is to set the Criteria.
  /// Configuration and Potentials (or System) must be set first.
  void set(std::shared_ptr<Criteria> criteria);

  // HWH depreciate
  void add(std::shared_ptr<Criteria> criteria) { set(criteria); }

  /// Once Criteria is set, it may be accessed on a read-only basis.
  const Criteria& criteria() const { return const_cast<Criteria&>(*criteria_); }

  /// Initialize the current energy.
  void initialize_energy();

  /// The remaining actions can be done in almost any order.
  /// Typically, one begins by adding trials.
  void add(std::shared_ptr<Trial> trial);

  /// Access the trials on a read-only basis.
  const TrialFactory& trials() const { return trial_factory_; }

  /// Access the trials on a read-only basis.
  const Trial& trial(const int index) const {
    return trial_factory_.trial(index); }

//  // HWH consider moving this to utils
//  // HWH consider depreciating this. While convenient, it confuses many users.
//  /// Before the production simulation, one may wish to quickly seek an initial
//  /// configuration with a given number of particles, without necessarily
//  /// satisfying detailed balance.
//  /// Thus, all stats are reset.
//  void seek_num_particles(const int num);

  /// An Analyzer performs some task after a given number of steps, but is
  /// read-only on System, Criteria and Trials.
  /// At this stage, multistate non-factory classes are converted into
  /// factories for each state in criteria.
  void add(std::shared_ptr<Analyze> analyze);

  /// Return all analyzers.
  const std::vector<std::shared_ptr<Analyze> >& analyzers() const {
    return analyze_factory_.analyzers(); }

  /// Return an Analyze by index.
  const Analyze& analyze(const int index) const {
    return analyze_factory_.analyze(index); }

  /// Return the number of analyzers.
  int num_analyzers() const {
    return static_cast<int>(analyze_factory_.analyzers().size()); }

  /// A Modifier performs some task after a given number of steps, but may
  /// change the System, Criteria and Trials.
  void add(const std::shared_ptr<Modify> modify);

  /// Return an Modify by index.
  const Modify& modify(const int index) const {
    return modify_factory_.modify(index); }

  /// Return the number of modifiers.
  int num_modifiers() const {
    return static_cast<int>(modify_factory_.modifiers().size()); }

  /// Add a checkpoint.
  void add(const std::shared_ptr<Checkpoint> checkpoint);

  /// Attempt one trial, with subsequent analysers and modifiers.
  // void attempt() { attempt_(1, &trial_factory_, random_.get()); }

  /// Attempt a number of Monte Carlo trials.
  void attempt(const int num_trials = 1) {
    attempt_(num_trials, &trial_factory_, random_.get()); }

  /// Reset trial statistics
  virtual void reset_trial_stats() { trial_factory_.reset_stats(); }

  /// Attempt Monte Carlo trials until Criteria returns completion.
  void run_until_complete() {
    run_until_complete_(&trial_factory_, random_.get()); }

  // HWH hackish interface for prefetch
  void before_attempts_();
  void delay_finalize_() {
    trial_factory_.delay_finalize(); }
  void after_trial_analyze_() {
    analyze_factory_.trial(*criteria_, system_, trial_factory_); }
  void after_trial_modify_();
  // Mimic a rejection by a trial.
  void imitate_trial_rejection_(const int trial_index,
      const double ln_prob,
      const int state_old,
      const int state_new);
  /// Attempt trial index without analyzers, modifiers or checkpoints.
  bool attempt_trial(const int index);
  // Revert changes from previous trial.
  void revert_(const int trial_index,
      const bool accepted,
      const double ln_prob);
  // Finalize changes from previous trial.
  void finalize_(const int trial_index) {
    trial_factory_.finalize(trial_index, &system_); }
  // Load random numbers and energy calculations into cache.
  void load_cache_(const bool load);
  // Unload random numbers and energy calculations from cache.
  void unload_cache_(const MonteCarlo& mc);

  virtual void serialize(std::ostream& ostr) const;
  explicit MonteCarlo(std::istream& istr);

  // HWH python interface cannot handle stringstreams with serialization.
  std::string serialize() {
    std::stringstream ss;
    serialize(ss);
    return ss.str();
  }
  MonteCarlo deserialize(const std::string str) {
    std::stringstream ss(str);
    return MonteCarlo(ss);
  }

  virtual ~MonteCarlo() {}

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

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_MONTE_CARLO_H_

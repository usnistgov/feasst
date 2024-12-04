
#ifndef FEASST_MONTE_CARLO_MONTE_CARLO_H_
#define FEASST_MONTE_CARLO_MONTE_CARLO_H_

#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <map>

namespace feasst {

class Action;
class Analyze;
class AnalyzeFactory;
class Checkpoint;
class Configuration;
class Criteria;
class Modify;
class ModifyFactory;
class NeighborCriteria;
class Potential;
class Random;
class Select;
class System;
class TimerRDTSC;
class ThermoParams;
class Trial;
class TrialFactory;
class TrialFactoryNamed;

typedef std::map<std::string, std::string> argtype;
typedef std::vector<std::pair<std::string, argtype> > arglist;

// HWH consider a constructor-based initialization of MonteCarlo..
// HWH something where order doesn't need to be enforced?
/**
  MonteCarlo contains Trials which perturb the System by generating the
  probability of acceptance through Criteria that are accepted based on a
  Random number generator.
  Between trials, MonteCarlo also contains classes that Analyze or Modify.
 */
class MonteCarlo {
 public:
  /// Construct with Random number generator.
  explicit MonteCarlo(std::shared_ptr<Random> random);

  /// Construct with RandomMT19937.
  MonteCarlo();

  /*
    Objects are processed with the following (derived-)names and arguments.
    - Random (default: RandomMT19937).
    - Configuration (multiple)
    - Potential (multiple; with configuration_index)
    - ThermoParams
    - Criteria
    - Trial (multiple)
    - Analyze (multiple)
    - Modify (multiple)
    - Action
    Those not labeled multiple will override any previous object of same type,
    while those labeled multiple will not override but rather add more.

    For example, in C++:

    auto mc = MakeMonteCarlo({{
      {"RandomMT19937", {{"seed", "123"}}},
      ...
    }});

    Or in Python:

    mc = fst.MakeMonteCarlo(fst.arglist([[
      "RandomMT19937", {"seed": "123"},
      ...
    ]]));

   */
  explicit MonteCarlo(arglist args, const bool silent = false);

  /// Begin processing the arguments.
  void begin(arglist args, const bool silent = false);

  /// Process arguments.
  void parse_args(arglist * args, const bool silent = false);

  /// Add more arguments.
  void add_args(arglist args);

  /// Resume processing the above arguments after Checkpointing.
  void resume(const bool silent = false);

  /// Clear arguments;
  void clear_arguments() { action_ = NULL; args_.clear(); }

  /// Set the random number generator.
  void set(std::shared_ptr<Random> random) { random_ = random; }

  /// Return the random number generator.
  const Random& random() const { return const_cast<Random&>(*random_); }

  /// Seed random number generator.
  void seed_random(const int seed);

  /// The first action with a Monte Carlo object is to set the Configuration.
  void add(std::shared_ptr<Configuration> config);

  /// The configuration may be accessed read-only.
  const Configuration& configuration(const int index = 0) const;

  /// The second action is to add Potentials.
  void add(std::shared_ptr<Potential> potential, const int config = 0);

  /// Warning for deprecated use.
  void add(const Potential& potential);

  /// Set an existing potential.
  void set(const int index, std::shared_ptr<Potential> potential);

  /// Add potential to optimized.
  void add_to_optimized(std::shared_ptr<Potential> potential);

  /// Add potential to reference.
  void add_to_reference(std::shared_ptr<Potential> potential,
    /// Store different references by index.
    const int index = 0,
    const int config = 0);

  /// Add NeighborCriteria.
  void add(std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const int config = 0);

  /// The third action is to set the ThermoParams.
  void set(std::shared_ptr<ThermoParams> thermo_params);

  /// Return the ThermoParams.
  const ThermoParams& thermo_params() const;

  /// Alternatively, the first, second and third actions may be combined by
  /// setting the system directly.
  /// This must be done before setting Criteria.
  void set(const System& system);

  /// Once the System is set, it may be accessed on a read-only basis.
  const System& system() const;

  /// Reinitialize the system. Return total energy.
  double initialize_system(const int config);

  /// Return the trials
  const TrialFactory& trial_factory() const { return *trial_factory_; }

  /// Return the analyzers
  const AnalyzeFactory& analyze_factory() const { return *analyze_factory_; }

  /// Return the modifiers
  const ModifyFactory& modify_factory() const { return *modify_factory_; }

  // HWH depreciate: only in rare cases should the system be modified directly.
  System * get_system();
  Criteria * get_criteria();
  Random * get_random() { return random_.get(); }
  TrialFactory * get_trial_factory();
  AnalyzeFactory * get_analyze_factory();
  ModifyFactory * get_modify_factory();

  // HWH hackish interface. See CollectionMatrixSplice::adjust_bounds.
  void adjust_bounds(const bool left_most, const bool right_most,
    const bool left_complete, const bool right_complete,
    const bool all_min_size,
    const int min_size, MonteCarlo * mc);

  /// The fourth action is to set the Criteria.
  /// Configuration and Potentials (or System) must be set first.
  void set(std::shared_ptr<Criteria> criteria);

  /// Once Criteria is set, it may be accessed on a read-only basis.
  const Criteria& criteria() const;

  /// Initialize the criteria. Also initializes system.
  void initialize_criteria();

  /// The remaining actions can be done in almost any order.
  /// Typically, one begins by adding trials.
  void add(std::shared_ptr<Trial> trial);

  /// Some Trials are simply TrialFactories containing multiple trials, such as
  /// TrialTransfer (factory of TrialAdd and TrialRemove).
  /// If a TrialFactory is added, flatten by adding individual trials instead.
  /// This means that add(MakeTrialTransfer()) will add both a TrialAdd and
  /// TrialRemove with weights equal to TrialFactory weight divided according to
  /// the weights of the individual trials.
  /// Thus, adding TrialTransfer with a weight of 4 will result in TrialAdd
  /// with weight of 2 and TrialRemove with weight of 2.
  void add(std::shared_ptr<TrialFactoryNamed> trials);

  /// Remove a trial by index.
  void remove_trial(const int index);

  /// Access the trials on a read-only basis.
  const TrialFactory& trials() const;

  /// Access the trials on a read-only basis.
  const Trial& trial(const int index) const;

  /// Initialize trials.
  void initialize_trials();

  /**
    An Analyzer performs some task after a given number of steps, but is
    read-only on System, Criteria and Trials.
    At this stage, multistate non-factory classes are converted into
    factories for each state in criteria.
   */
  void add(std::shared_ptr<Analyze> analyze);

  /// Remove an analyze by index.
  void remove_analyze(const int index);

  /// Return all analyzers.
  const std::vector<std::shared_ptr<Analyze> >& analyzers() const;

  /// Return an Analyze by index.
  const Analyze& analyze(const int index) const;

  /// Return the number of analyzers.
  int num_analyzers() const;

  /// Initialize analyzers.
  void initialize_analyzers();

  /// A Modifier performs some task after a given number of steps, but may
  /// change the System, Criteria and Trials.
  void add(const std::shared_ptr<Modify> modify);

  /// Remove a modify by index.
  void remove_modify(const int index);

  /// Return all modifiers.
  const std::vector<std::shared_ptr<Modify> >& modifiers() const;

  /// Return an Modify by index.
  const Modify& modify(const int index) const;

  /// Return the number of modifiers.
  int num_modifiers() const;

  /// Add a checkpoint.
  void set(const std::shared_ptr<Checkpoint> checkpoint);

  /// Write checkpoint file
  void write_checkpoint() const;

  /// Attempt one trial, with subsequent analysers and modifiers.
  // void attempt() { attempt_(1, trial_factory_.get(), random_.get()); }

  /// Perform an Action
  virtual void run(std::shared_ptr<Action> action);

  /// Attempt a number of Monte Carlo trials.
  void attempt(const int num_trials = 1);

  /// Reset trial statistics
  virtual void reset_trial_stats();

  /// Run a number of trials.
  virtual void run_num_trials(int num_trials);

  /// Run until a number of particles is reached.
  virtual void run_until_num_particles(const int num_particles,
                                       const int particle_type,
                                       const int configuration_index);

  /// Run trials for a number of hours.
  virtual void run_for_hours(const double hours);

  /// Attempt Monte Carlo trials until Criteria returns completion.
  /// If available, automatically write checkpoint when complete.
  void run_until_complete();

  /// Attempt Monte Carlo trials until the given file name exists.
  virtual void run_until_file_exists(const std::string& file_name);

  /// Set the Criteria::cycles_to_complete.
  void set_cycles_to_complete(const int num);

  // HWH hackish interface for prefetch
  void before_attempts_();
  void delay_finalize_();
  void after_trial_analyze_();
  void after_trial_modify_();
  void after_trial_checkpoint_();
  // Mimic a rejection by a trial.
  void imitate_trial_rejection_(const int trial_index,
    const double ln_prob,
    const bool endpoint,
    const bool auto_reject,
    const int state_old,
    const int state_new);
  void ghost_trial_(
    const double ln_prob,
    const int state_old,
    const int state_new,
    const bool endpoint);
  /// Attempt trial index without analyzers, modifiers or checkpoints.
  bool attempt_trial(const int index);
  // Revert changes from previous trial.
  void revert_(const int trial_index,
    const bool accepted,
    const bool endpoint,
    const bool auto_reject,
    const double ln_prob);
  // Finalize changes from previous trial.
  void finalize_(const int trial_index);
  // Load random numbers and energy calculations into cache.
  void load_cache_(const bool load);
  // Unload random numbers and energy calculations from cache.
  void unload_cache_(const MonteCarlo& mc);
  void synchronize_(const MonteCarlo& mc, const Select& perturbed);

  /// Write all Analyze and Modify.
  void write_to_file();

  /// Return the last argument parsed.
  const std::pair<std::string, argtype>& next_arg() const {
    return next_arg_; }

  /// Set the number of configurations to parse simultaneously
  void set_parse_for_num_configs(const int num);

  /// Set the parse replacements
  void set_parse_replace(
    const std::vector<std::vector<std::string> >& replace = {});

  /// Set the parse replacement
  void set_replace_with_index(const std::string& str);

  /// Set the timer
  void set_timer();

  /// Return timer
  const TimerRDTSC * const timer() const { return timer_.get(); }

  virtual void serialize(std::ostream& ostr) const;
  explicit MonteCarlo(std::istream& istr);

  // HWH python interface cannot handle stringstreams with serialization.
  std::string serialize() const;
//  MonteCarlo deserialize(const std::string str);

  virtual ~MonteCarlo();

 protected:
  virtual void attempt_(int num_trials,
    TrialFactory * trial_factory,
    Random * random);
  virtual void run_until_complete_(TrialFactory * trial_factory,
                                   Random * random);

 private:
  std::unique_ptr<System> system_;
  std::shared_ptr<Criteria> criteria_;
  std::unique_ptr<TrialFactory> trial_factory_;
  std::unique_ptr<AnalyzeFactory> analyze_factory_;
  std::unique_ptr<ModifyFactory> modify_factory_;
  std::shared_ptr<Checkpoint> checkpoint_;
  std::shared_ptr<Random> random_;
  std::shared_ptr<Action> action_;
  arglist args_;

  bool config_set_ = false;
  bool potential_set_ = false;
  bool thermo_params_set_ = false;
  bool system_set_ = false;
  bool criteria_set_ = false;
  int parse_for_num_configs_ = 1;
  std::vector<std::vector<std::string> > parse_replace_;
  std::string replace_with_index_;

  std::unique_ptr<TimerRDTSC> timer_;

  // temporary and not serialized
  std::pair<std::string, argtype> next_arg_;

  void record_next_arg_(arglist * args);
  bool duplicate_stepper_output_file_(const std::string output_file);
};

inline std::shared_ptr<MonteCarlo> MakeMonteCarlo() {
  return std::make_shared<MonteCarlo>();
}

inline std::shared_ptr<MonteCarlo> MakeMonteCarlo(arglist args,
    const bool silent = false) {
  return std::make_shared<MonteCarlo>(args, silent); }

/// Construct MonteCarlo from file.
std::shared_ptr<MonteCarlo> MakeMonteCarlo(const std::string file_name);

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_MONTE_CARLO_H_

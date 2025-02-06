
#ifndef FEASST_MONTE_CARLO_CRITERIA_H_
#define FEASST_MONTE_CARLO_CRITERIA_H_

#include <vector>
#include <map>
#include <string>
#include <memory>
#include "system/include/synchronize_data.h"

namespace feasst {

class Acceptance;
class Bias;
class BondVisitor;
class Constraint;
class FlatHistogram;
class Macrostate;
class Random;
class System;

typedef std::map<std::string, std::string> argtype;

/**
  Determine whether to accept or reject a trial.
  Stores the total energy based on energy changes from each trial.
 */
class Criteria {
 public:
  //@{
  /** @name Arguments
    - cycles_to_complete: set the number of cycles for a simulation
      to be considered complete (default: 20).
    - Constraint: ConstrainNumParticles, AHalfB, etc.
   */
  explicit Criteria(argtype args = argtype());
  explicit Criteria(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Same as above, but also add a constraint.
  explicit Criteria(std::shared_ptr<Constraint> constraint,
                    argtype args = argtype());

  /// Add a constraint.
  void add(std::shared_ptr<Constraint> constraint) {
    constraints_.push_back(constraint); }

  virtual void precompute(System * system) {}

  /// Return whether constraints are statisfied.
  bool is_allowed(const System& system,
                  const Acceptance& acceptance);

  /// This function is called before a trial attempt.
  virtual void before_attempt(const System& system) {}

  /// Return whether or not the trial attempt should be accepted.
  virtual bool is_accepted(
    const System& system,
    Acceptance * acceptance,
    Random * random) = 0;

  virtual void finalize(const Acceptance& acceptance) {}
  virtual void revert(const Acceptance& acceptance) {}

  /// Return whether or not the last trial attempt was accepted.
  bool was_accepted() const { return was_accepted_; }
  void set_was_accepted(const bool acc) { was_accepted_ = acc; }

  /// Set the current total energy based on energy changes per trial in order
  /// to avoid recomputation of the energy of the entire configuration.
  /// For example, Metropolis Monte Carlo trials are concerned with the change
  /// in energy, and this variable tracks the total from the changes.
  void set_current_energy(const double energy, const int config = 0);

  /// Return the current total energy based on energy changes per trial.
  double current_energy(const int config = 0) const;

  /// Same as above, except for energy profiles instead of total energy.
  void set_current_energy_profile(const std::vector<double>& energy,
    const int config = 0);

  /// Return the current energy profile based on energy changes per trial.
  const std::vector<double>& current_energy_profile(const int config = 0) const;

  /// Update the current energy.
  void update_current_energy(const Acceptance& acceptance);

  /// Return the header of the status for periodic output.
  std::string status_header(const System& system,
                            const bool include_bonds) const;

  /// Return the brief status for periodic output.
  std::string status(const System& system, const bool max_precision,
    const bool include_bonds, BondVisitor * visitor) const;

  /// Return a human-readable output of all data (not as brief as status).
  virtual std::string write() const;

  /// Return the simulation phase index used to differentiate production
  /// and initialization, etc.
  virtual int phase() const { return phase_; }

  /// Increment the simulation phase.
  virtual void increment_phase() { ++phase_; }

  /// Return the number of cycles for a simulation to be complete.
  /// Iterations are defined by the Derived class.
  /// For example, one Metropolis cycle is 1000 trials.
  /// FlatHistogram cycles depend on the Bias.
  /// For TransitionMatrix, one cycle is a sweep.
  virtual int cycles_to_complete() const {
    return cycles_to_complete_; }

  /// Set the number of cycles for a simulation to be complete.
  virtual void set_cycles_to_complete(const int num) {
    cycles_to_complete_ = num; }

  /// Return the current number of cycles.
  virtual int num_cycles(
    /// If != -1, return cycles of a particular state (TM/WLTM only).
    const int state = -1) const;

  /// Return true if the number of cycles for completion has been reached.
  virtual bool is_complete() const {
    return num_cycles() >= cycles_to_complete(); }

  /// Set the simulation as complete. Used for post processing.
  virtual void set_complete() { set_cycles_to_complete(0); }

  /// Return the state index for multistate simulations (default: 0).
  virtual int state() const { return 0; }

  /// Return the number of states. (default: 1).
  virtual int num_states() const { return 1; }

  // HWH deprecate this. Only used by old growth expanded, not morph.
  /// Set the expanded state.
  void set_expanded_state(const int state = 0, const int num = 1);

  /// Return the expanded state.
  int expanded_state() const { return expanded_state_; }

  /// Return number of expanded states.
  int num_expanded_states() const { return num_expanded_states_; }

  // HWH for prefetch
  virtual int state_old() const { return 0; }
  virtual int state_new() const { return 0; }

  /// Update.
  virtual void update() {}

  /// Return true if equivalent.
  virtual bool is_equal(const Criteria& criteria,
                        const double tolerance) const;

  /// Same as above but with a default tolerance of NEAR_ZERO
  virtual bool is_equal(const Criteria& criteria) const;

  // HWH hackish interface for prefetch
  // Revert changes from previous trial.
  virtual void revert_(const bool accepted, const bool endpoint,
                       const double ln_prob, const std::vector<int>& updated);
  // Imitate a trial rejection (used in FlatHistogram).
  virtual void imitate_trial_rejection_(const double ln_prob,
    const int state_old,
    const int state_new,
    const bool endpoint) {}
  void synchronize_(const Criteria& criteria);
  const SynchronizeData& data() const { return data_; }

  // HWH hackish adjust_bounds interface. See CollectionMatrixSplice.
  virtual int set_soft_max(const int index, const System& sys);
  virtual int set_soft_min(const int index, const System& sys);
  virtual void set_cm(const bool inc_max, const int macro,
                      const Criteria& crit);
  virtual void adjust_bounds(const bool left_most, const bool right_most,
    const bool left_complete, const bool right_complete,
    const bool all_min_size,
    const int min_size, const System& system, const System * upper_sys,
    Criteria * criteria, bool * adjusted_up, std::vector<int> * states);
  virtual const Macrostate& macrostate() const;
  virtual int soft_max() const;
  virtual int soft_min() const;
  virtual const Bias& bias() const;
  virtual const FlatHistogram& flat_histogram() const;

  // HWH hackish interface for setting the state in post processing.
  virtual void update_state(const System& system, const Acceptance& accept) {}

  /// Initialize criteria
  void initialize(System * system);

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Criteria> create(std::istream& istr) const;
  virtual std::shared_ptr<Criteria> create(argtype * args) const;
  std::map<std::string, std::shared_ptr<Criteria> >& deserialize_map();
  std::shared_ptr<Criteria> deserialize(std::istream& istr);
  std::shared_ptr<Criteria> factory(const std::string name, argtype * args);
  explicit Criteria(std::istream& istr);
  virtual ~Criteria() {}

  //@}
 protected:
  std::string class_name_ = "Criteria";
  void serialize_criteria_(std::ostream& ostr) const;
  bool was_accepted_ = false;
  SynchronizeData data_;
  void check_num_cycles_(const int trials_per_cycle);

 private:
  std::vector<double> previous_energy_;
  std::vector<std::vector<double> > previous_energy_profile_;
  int phase_ = 0;
  int expanded_state_;
  int num_expanded_states_;
  int cycles_to_complete_;
  std::vector<std::shared_ptr<Constraint> > constraints_;

  std::vector<double> * current_energy_() {
    return &((*data_.get_dble_2D())[0]); }
  const std::vector<double>& const_current_energy_() const {
    return data_.dble_2D()[0]; }
  std::vector<std::vector<double> > * current_energy_profile_() {
    return &((*data_.get_dble_3D())[0]); }
  const std::vector<std::vector<double> >& const_current_energy_profile_()
    const { return data_.dble_3D()[0]; }
  int * num_attempt_since_last_cycle_() {
    return &((*data_.get_int_1D())[0]); }
  int * num_cycles_() { return &((*data_.get_int_1D())[1]); }
  int const_num_cycles_() const { return data_.int_1D()[1]; }
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_CRITERIA_H_

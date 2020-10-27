
#ifndef FEASST_MONTE_CARLO_CRITERIA_H_
#define FEASST_MONTE_CARLO_CRITERIA_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "system/include/synchronize_data.h"
#include "monte_carlo/include/acceptance.h"

namespace feasst {

class Random;
class Constraint;

/**
  Determine whether to accept or reject a trial.
  Stores the total energy based on energy changes from each trial.
  probability of acceptance.
 */
class Criteria {
 public:
  explicit Criteria(const argtype& args = argtype());

  /// Same as above, but also add a constraint.
  Criteria(std::shared_ptr<Constraint> constraint,
    const argtype& args = argtype());

  /// Add a constraint.
  void add(std::shared_ptr<Constraint> constraint) {
    constraints_.push_back(constraint); }

  /// Return whether constraints are statisfied.
  bool is_allowed(const System& system,
                  const Acceptance& acceptance);

  /// This function is called before a trial attempt.
  virtual void before_attempt(const System& system) {}

  /// Return whether or not the trial attempt should be accepted.
  virtual bool is_accepted(const Acceptance& acceptance_,
    const System& system,
    Random * random) = 0;

  /// Return whether or not the last trial attempt was accepted.
  // bool was_accepted() const { return was_accepted_; }

  /// Set the current total energy based on energy changes per trial in order
  /// to avoid recomputation of the energy of the entire configuration.
  /// For example, Metropolis Monte Carlo trials are concerned with the change
  /// in energy, and this variable tracks the total from the changes.
  void set_current_energy(const double energy);

  /// Return the current total energy based on energy changes per trial.
  double current_energy() const { return data_.dble_1D()[0]; }

  /// Return the header of the status for periodic output.
  std::string status_header() const;

  /// Return the brief status for periodic output.
  std::string status() const;

  /// Return a human-readable output of all data (not as brief as status).
  virtual std::string write() const;

  /// Return true if completion requirements are met.
  virtual bool is_complete() const { return false; }

  /// Return the simulation phase index used to differentiate production
  /// and initialization, etc.
  virtual int phase() const { return phase_; }

  /// Increment the simulation phase.
  virtual void increment_phase() { ++phase_; }

  // HWH consider using this to set number of trials in Metropolis
  /// Set the number of iterations.
  virtual void set_num_iterations(const int iteration);

  /// Return the state index for multistate simulations (default: 0).
  virtual int state() const { return 0; }

  /// Return the number of states. (default: 1).
  virtual int num_states() const { return 1; }

  // HWH depreciate this. Only used by old growth expanded, not morph.
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
  virtual void revert_(const bool accepted, const double ln_prob);
  // Imitate a trial rejection (used in FlatHistogram).
  virtual void imitate_trial_rejection_(const double ln_prob,
    const int state_old,
    const int state_new) {}
  void synchronize_(const Criteria& criteria) { data_ = criteria.data(); }
  const SynchronizeData& data() const { return data_; }

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Criteria> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Criteria> >& deserialize_map();
  std::shared_ptr<Criteria> deserialize(std::istream& istr);
  Criteria(std::istream& istr);
  virtual ~Criteria() {}

 protected:
  std::string class_name_ = "Criteria";
  Arguments args_;
  void serialize_criteria_(std::ostream& ostr) const;
  bool was_accepted_ = false;
  SynchronizeData data_;

 private:
  double * current_energy_() { return &((*data_.get_dble_1D())[0]); }
  double previous_energy_ = 0.;
  int phase_ = 0;
  int expanded_state_;
  int num_expanded_states_;
  std::vector<std::shared_ptr<Constraint> > constraints_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_CRITERIA_H_

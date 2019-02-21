
#ifndef FEASST_CORE_CRITERIA_H_
#define FEASST_CORE_CRITERIA_H_

#include <vector>
#include "core/include/system.h"

namespace feasst {

// This is the information passed from the system to the acceptance criteria in
// order to determine whether to accept or reject the trial.
struct AcceptanceCriteria {
  double ln_metropolis_prob = 0.;
  double energy_new = 0.;
  double energy_new_select = 0.;
  int force_rejection = 0;
  System* system = NULL;
  int accepted = -1;
};

/**
  Determine whether to accept or reject a trial.
  Stores the total energy based on energy changes from each trial.
  Stores thermodynamic variables relevant to computing the Metropolis
  probability of acceptance.
 */
class Criteria {
 public:
  /// Set beta, the inverse temperature \f$ \beta=\frac{1}{k_B T} \f$.
  void set_beta(const double beta);

  /// Return beta.
  double beta() const;

  /// Add an activity for a given type of particle.
  /// Note that z has units length^{-dimension} such that Vz/N is unitless.
  // HWH Note: consider interfacing somehow with system/particle, etc
  // Perhaps System should contain criteria or MC kernel new object
  void add_activity(const double activity) { activity_.push_back(activity); }

  /// Return the activity of the particle type.
  double activity(const int particle_type = 0) const;

  /// This function is called before a trial attempt.
  virtual void before_attempt(const System* system) {}

  /// Return whether or not the trial attempt should be accepted.
  virtual bool is_accepted(const AcceptanceCriteria accept_criteria) = 0;

  /// Set the current total energy based on energy changes per trial in order
  /// to avoid avoid recomputation of the energy of the entire configuration.
  /// For example, Metropolis Monte Carlo trials are concerned with the change
  /// in energy, and this variable tracks the total from the changes.
  void set_running_energy(const double energy) { running_energy_ = energy; }

  /// Return the current total energy based on energy changes per trial.
  double running_energy() { return running_energy_; }

  virtual ~Criteria() {}

 private:
  double beta_;
  bool beta_initialized_ = false;
  std::vector<double> activity_;
  double running_energy_;

  /// This function is called after a trial attempt but before acceptance
  /// decision.
  virtual void after_attempt_(const System* system) {}
};

}  // namespace feasst

#endif  // FEASST_CORE_CRITERIA_H_


#ifndef FEASST_MONTE_CARLO_CRITERIA_H_
#define FEASST_MONTE_CARLO_CRITERIA_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/acceptance.h"

namespace feasst {

/**
  Determine whether to accept or reject a trial.
  Stores the total energy based on energy changes from each trial.
  Stores thermodynamic variables relevant to computing the Metropolis
  probability of acceptance.
 */
class Criteria {
 public:
  Criteria(
    /**
      beta : inverse temperature, \f$ \beta = \frac{1}{k_B T} \f$.

      chemical_potential[i] : chemical potential of the i-th particle type.
        The [i] is to be substituted for an integer 0, 1, 2, ...
        If only one particle type, you can drop the [i].
        The chemical potential must have the inverse units of \f$\beta\f$.
     */
    const argtype& args = argtype());

  /// Set beta, the inverse temperature \f$ \beta=\frac{1}{k_B T} \f$.
  void set_beta(const double beta);

  /// Return beta.
  double beta() const;

  /// Add a chemical potential for a given type of particle.
  /// Note that z has units length^{-dimension} such that Vz/N is unitless.
  // HWH Note: consider interfacing somehow with system/particle, etc
  // Perhaps System should contain criteria or MC kernel new object
  void add_chemical_potential(const double chemical_potential) {
    chemical_potentials_.push_back(chemical_potential); }

  /// Set the chemical potential of a given type.
  void set_chemical_potential(const double mu, const int particle_type = 0) {
    chemical_potentials_[particle_type] = mu; }

  /// Return the chemical potential of the particle type.
  double chemical_potential(const int particle_type = 0) const;

  /// Return the dimensionless product of beta and the chemical potential.
  double beta_mu(const int particle_type = 0) const {
    return beta()*chemical_potentials_[particle_type]; }

  /// This function is called before a trial attempt.
  virtual void before_attempt(const System* system) {}

  /// Return whether or not the trial attempt should be accepted.
  virtual bool is_accepted(const Acceptance& acceptance_,
    const System * system,
    const double uniform_random) = 0;

  /// Set the current total energy based on energy changes per trial in order
  /// to avoid recomputation of the energy of the entire configuration.
  /// For example, Metropolis Monte Carlo trials are concerned with the change
  /// in energy, and this variable tracks the total from the changes.
  void set_current_energy(const double energy) {
    previous_energy_ = current_energy_;
    current_energy_ = energy;
  }

  /// Return the current total energy based on energy changes per trial.
  double current_energy() const { return current_energy_; }

  /// Return the header of the status for periodic output.
  std::string status_header() const;

  /// Return the brief status for periodic output.
  std::string status() const;

  /// Return a human-readable output of all data (not as brief as status).
  virtual std::string write() const;

  /// Return true if completion requirements are met.
  virtual bool is_complete() { return false; }

  /// Return the state index for multistate simulations (default: 0).
  virtual int state() const { return 0; }

  /// Return the number of states. (default: 1).
  virtual int num_states() const { return 1; }

  /// Set the trial state.
  void set_trial_state(const int state = 0, const int num = 1);

  /// Return the trial state.
  int trial_state() const { return trial_state_; }

  /// Return number of trial states.
  int num_trial_states() const { return num_trial_states_; }

  /// Revert changes from previous trial.
  virtual void revert(const bool accepted);

  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Criteria> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Criteria> >& deserialize_map();
  std::shared_ptr<Criteria> deserialize(std::istream& istr);
  Criteria(std::istream& istr);
  virtual ~Criteria() {}

 protected:
  Arguments args_;
  void serialize_criteria_(std::ostream& ostr) const;

 private:
  double beta_;
  bool beta_initialized_ = false;
  std::vector<double> chemical_potentials_;
  double current_energy_ = 0.;
  double previous_energy_ = 0.;
  int trial_state_;
  int num_trial_states_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_CRITERIA_H_

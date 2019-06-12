
#ifndef FEASST_MONTE_CARLO_CRITERIA_H_
#define FEASST_MONTE_CARLO_CRITERIA_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"

namespace feasst {

/*
  This object contains information necessary for Criteria to make a decision on
  whether or not to accept or reject a trial.
 */
class Acceptance {
 public:
  Acceptance() { reset(); }

  /// Return the natural logarithm of the Metropolis acceptance probability.
  double ln_metropolis_prob() const { return ln_metropolis_prob_; }

  /// Set the above quantity.
  void set_ln_metropolis_prob(const double prob = 0.) {
    ln_metropolis_prob_ = prob; }

  /// Add to the above quantity.
  void add_to_ln_metropolis_prob(const double prob = 0.) {
    ln_metropolis_prob_ += prob; }

  /// Return whether or not to reject the trial outright.
  bool reject() const { return reject_; }

  /// Set the above quantity.
  void set_reject(const bool reject = false) {
    reject_ = reject; }

  /// Reset all stored quantities before each trial.
  void reset() {
    set_ln_metropolis_prob();
    set_reject();
    energy_new_ = 0.;
    energy_old_ = 0.;
    macrostate_shift_ = 0;
  }

  /// Return the energy of the new configuration.
  double energy_new() const { return energy_new_; }

  /// Set the above quantity.
  void set_energy_new(const double energy) { energy_new_ = energy; }

  /// Add to the above quantity.
  void add_to_energy_new(const double energy) { energy_new_ += energy; }

  /// Return the energy of the old configuration.
  double energy_old() const { return energy_old_; }

  /// Set the above quantity.
  void set_energy_old(const double energy) { energy_old_ = energy; }

  /// Return the above quantity.
  void add_to_energy_old(const double energy) { energy_old_ += energy; }

  /// Return the energy of the reference.
  double energy_ref() const { return energy_ref_; }

  /// Set the above quantity.
  void set_energy_ref(const double energy) { energy_ref_ = energy; }

  /// Return the shift in the macrostate due to an optimization where
  /// Perturb does not completely update system until finalize.
  int macrostate_shift() const { return macrostate_shift_; }

  /// Set the above.
  void set_macrostate_shift(const int shift) { macrostate_shift_ = shift; }

 private:
  double ln_metropolis_prob_;
  double energy_new_;
  double energy_old_;
  double energy_ref_;
  int macrostate_shift_;
  bool reject_;
};

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

  /// Return the chemical potential of the particle type.
  double chemical_potential(const int particle_type = 0) const;

  /// Return the dimensionless product of beta and the chemical potential.
  double beta_mu(const int particle_type = 0) const {
    return beta()*chemical_potentials_[particle_type]; }

  /// This function is called before a trial attempt.
  virtual void before_attempt(const System* system) {}

  /// Return whether or not the trial attempt should be accepted.
  virtual bool is_accepted(const Acceptance& acceptance_,
    const System * system) = 0;

  /// Set the current total energy based on energy changes per trial in order
  /// to avoid recomputation of the energy of the entire configuration.
  /// For example, Metropolis Monte Carlo trials are concerned with the change
  /// in energy, and this variable tracks the total from the changes.
  void set_current_energy(const double energy) { current_energy_ = energy; }

  /// Return the current total energy based on energy changes per trial.
  double current_energy() const { return current_energy_; }

  /// Return the header of the status for periodic output.
  std::string status_header() const {
    return std::string("energy"); }

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
  double current_energy_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_CRITERIA_H_

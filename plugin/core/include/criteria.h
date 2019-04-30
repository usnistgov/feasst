
#ifndef FEASST_CORE_CRITERIA_H_
#define FEASST_CORE_CRITERIA_H_

#include <vector>
#include "core/include/arguments.h"
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
  virtual bool is_accepted(const AcceptanceCriteria accept_criteria) = 0;

  /// Set the current total energy based on energy changes per trial in order
  /// to avoid recomputation of the energy of the entire configuration.
  /// For example, Metropolis Monte Carlo trials are concerned with the change
  /// in energy, and this variable tracks the total from the changes.
  void set_running_energy(const double energy) { running_energy_ = energy; }

  /// Return the current total energy based on energy changes per trial.
  double running_energy() const { return running_energy_; }

  /// Return the header of the status for periodic output.
  std::string status_header() const {
    return std::string("energy"); }

  /// Return the brief status for periodic output.
  std::string status() const;

  /// Return a human-readable output of all data (not as brief as status).
  virtual std::string write() const;

  /// Return true if completion requirements are met.
  virtual bool is_complete() { return false; }

  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Criteria> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Criteria> >& deserialize_map();
  std::shared_ptr<Criteria> deserialize(std::istream& istr);
  Criteria(std::istream& istr);
  virtual ~Criteria() {}

 protected:
  Arguments args_;
  void serialize_criteria_(std::ostream& ostr) const {
    feasst_serialize_version(1, ostr);
    feasst_serialize(beta_, ostr);
    feasst_serialize(beta_initialized_, ostr);
    feasst_serialize(chemical_potentials_, ostr);
    feasst_serialize(running_energy_, ostr);
  }

 private:
  double beta_;
  bool beta_initialized_ = false;
  std::vector<double> chemical_potentials_;
  double running_energy_;

  /// This function is called after a trial attempt but before acceptance
  /// decision.
  virtual void after_attempt_(const System* system) {}
};

}  // namespace feasst

#endif  // FEASST_CORE_CRITERIA_H_

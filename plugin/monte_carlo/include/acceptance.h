
#ifndef FEASST_MONTE_CARLO_ACCEPTANCE_H_
#define FEASST_MONTE_CARLO_ACCEPTANCE_H_

#include "configuration/include/select.h"

namespace feasst {

/**
  This object contains information necessary for Criteria to make a decision on
  whether or not to accept or reject a trial.
 */
class Acceptance {
 public:
  Acceptance() { reset(); }

  /// Return the natural logarithm of the Metropolis acceptance probability.
  double ln_metropolis_prob() const;

  /// Set the above quantity.
  void set_ln_metropolis_prob(const double prob = 0) {
    ln_metropolis_prob_ = prob; }

  /// Add to the above quantity.
  void add_to_ln_metropolis_prob(const double prob = 0) {
    ln_metropolis_prob_ += prob; }

  /// Return whether or not to reject the trial outright.
  bool reject() const { return reject_; }

  /// Set the above quantity.
  void set_reject(const bool reject = false) {
    reject_ = reject; }

  /// Reset all stored quantities before each trial.
  void reset();

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
  /// This assumes a particular macrostate is used for the given trial.
  // HWH refactor this so that perturb_remove shows deleted particle but
  // cell lists, etc are only updated upon finalization (optimization).
  int macrostate_shift() const { return macrostate_shift_; }
  int macrostate_shift_type() const { return macrostate_shift_type_; }

  /// Add to the above.
  void add_to_macrostate_shift(const int shift) { macrostate_shift_ += shift; }
  void set_macrostate_shift_type(const int type) {
    macrostate_shift_type_ += type; }

  /// Add to perturbed selection and equate trial state.
  void add_to_perturbed(const Select& select);

  /// Set perturbed trial state.
  void set_perturbed_state(const int state);

  /// Return the perturbed selection.
  const Select& perturbed() const { return perturbed_; }

 private:
  double ln_metropolis_prob_;
  double energy_new_;
  double energy_old_;
  double energy_ref_;
  int macrostate_shift_;
  int macrostate_shift_type_;
  bool reject_;
  Select perturbed_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_ACCEPTANCE_H_

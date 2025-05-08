
#ifndef FEASST_MONTE_CARLO_ACCEPTANCE_H_
#define FEASST_MONTE_CARLO_ACCEPTANCE_H_

#include <memory>
#include <vector>

namespace feasst {

class Select;

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
  void set_reject(const bool reject = false);

  /// Return whether or not to the trial was at a macrostate endpoint.
  bool endpoint() const { return endpoint_; }

  /// Set the above quantity.
  void set_endpoint(const bool endpoint = true) {
    endpoint_ = endpoint; }

  /// Reset all stored quantities before each trial.
  void reset();

  /// Return 1 if this conf has updated energy.
  int updated(const int conf = 0) const;
  const std::vector<int>& updtd() const { return updated_; }

  /// Return the energy of the new configuration.
  double energy_new(const int config = 0) const;

  /// Set the above quantity.
  void set_energy_new(const double energy, const int config = 0);

  /// Add to the above quantity.
  void add_to_energy_new(const double energy, const int config = 0);

  /// Return the configuration indices.
  int num_configurations() const;

  /// Return the energy profile of the new configuration.
  const std::vector<double>& energy_profile_new(const int config = 0) const;

  /// Set the above quantity.
  void set_energy_profile_new(const std::vector<double>& energy,
                              const int config = 0);

  /// Add to the above quantity.
  void add_to_energy_profile_new(const std::vector<double>& energy,
                                 const int config = 0);

  /// Subtract from the above quantity.
  void subtract_from_energy_profile_new(const std::vector<double>& energy,
                                        const int config = 0);

  /// Return the energy of the old configuration.
  double energy_old(const int config = 0) const;

  /// Set the above quantity.
  void set_energy_old(const double energy, const int config = 0);

  /// Add to the above quantity.
  void add_to_energy_old(const double energy, const int config = 0);

  /// Return the energy profile of the old configuration.
  const std::vector<double>& energy_profile_old(const int config = 0) const;

  /// Set the above quantity.
  void set_energy_profile_old(const std::vector<double>& energy,
                              const int config = 0);

  /// Add to the above quantity.
  void add_to_energy_profile_old(const std::vector<double>& energy,
                                 const int config = 0);

  /// Return the energy of the reference.
  double energy_ref(const int config = 0) const;

  /// Set the above quantity.
  void set_energy_ref(const double energy, const int config = 0);

  /// Return the shift in the macrostate due to an optimization where
  /// Perturb does not completely update system until finalize.
  /// This assumes a particular macrostate is used for the given trial.
  // HWH refactor this so that perturb_remove shows deleted particle but
  // cell lists, etc are only updated upon finalization (optimization).
  int macrostate_shift(const int config = 0) const {
    return macrostate_shift_[config]; }
  int macrostate_shift_type(const int config = 0) const {
    return macrostate_shift_type_[config]; }

  /// Add to the above.
  void add_to_macrostate_shift(const int shift, const int config = 0);
  void set_macrostate_shift_type(const int type, const int config = 0);

  /// Add to perturbed selection and equate trial state.
  void add_to_perturbed(const Select& select, const int config = 0);

  /// Set perturbed trial state.
  void set_perturbed_state(const int state, const int config = 0);

  void sort_perturbed();

  /// Return the perturbed selection.
  const std::vector<std::shared_ptr<Select> >& perturbed() const { return perturbed_; }
  const Select& perturbed(const int config) const;

  ~Acceptance();

 private:
  double ln_metropolis_prob_;
  std::vector<double> energy_new_;
  std::vector<double> energy_old_;
  std::vector<double> energy_ref_;
  std::vector<int> macrostate_shift_;
  std::vector<int> macrostate_shift_type_;
  bool reject_;
  bool endpoint_;
  std::vector<std::vector<double> > energy_profile_new_;
  std::vector<std::vector<double> > energy_profile_old_;
  std::vector<std::shared_ptr<Select> > perturbed_;
  std::vector<int> updated_;

  template <typename T>
  void resize_(const int config, std::vector<T> * vec) {
    if (config >= static_cast<int>(vec->size())) {
      vec->resize(config + 1);
    }
  }
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_ACCEPTANCE_H_

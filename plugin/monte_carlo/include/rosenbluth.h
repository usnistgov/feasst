
#ifndef FEASST_MONTE_CARLO_ROSENBLUTH_H_
#define FEASST_MONTE_CARLO_ROSENBLUTH_H_

#include "configuration/include/select.h"

namespace feasst {

/**
Store the energies and Rosenbluth factors necessary for configurational bias
Monte Carlo trials.

\f$w = \frac{1}{n} \sum_i^n \exp^{-\beta \Delta U}\f$

where \f$w\f$ is the Rosenbluth weight for \f$n\f$ steps.
 */
class Rosenbluth {
 public:
  Rosenbluth() {}

  /// Resize the energy, Boltmann factors by the number of steps in the trial.
  void resize(const int num);

  /// Return the number of steps.
  int num() const { return static_cast<int>(energy_.size()); }

  /// Store the selection for each step.
  void store(const int step, const Select& select) { stored_[step] = select; }

  /// Set the energy and excluded energy of the step.
  /// Exclude energy includes terms such as bond potentials where the bonds
  /// where already selected according to the appropriate distribution.
  void set_energy(const int step, const double energy, const double excluded);
  void set_energy_profile(const int step, const std::vector<double>& energy);

  /// Compute Boltzmann factors and cumulative probabilities for all steps.
  /// Choose one of the steps based on the probabilities.
  void compute(const double beta, Random * random, const bool old);

  /// Return the stored selection of the step.
  const Select& stored(const int step) const { return stored_[step]; }

  /// Return the chosen selection.
  const Select& chosen() const;

  /// Return the chosen energy.
  double chosen_energy() const;

  /// Return the chosen energy profile.
  const std::vector<double>& chosen_energy_profile() const;

  /// Return the natural logrithm of the total rosenbluth factor.
  double ln_total_rosenbluth() const { return ln_total_rosenbluth_; }

  /// Return the energy for a given step.
  double energy(const int step) const { return energy_[step]; }

  /// Return the energy profile for a given step.
  const std::vector<double> energy_profile(const int step) const {
    return energy_profile_[step]; }

  /// Return the chosen step
  int chosen_step() const { return chosen_step_; }

  /// Serialize.
  void serialize(std::ostream& ostr) const;

  /// Deserialize
  explicit Rosenbluth(std::istream& istr);

 private:
  std::vector<double> energy_;
  std::vector<std::vector<double> > energy_profile_;
  std::vector<double> excluded_;
  std::vector<double> weight_;
  std::vector<double> cumulative_;
  std::vector<Select> stored_;

  // temporary
  double ln_total_rosenbluth_;
  int chosen_step_ = -1;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_ROSENBLUTH_H_

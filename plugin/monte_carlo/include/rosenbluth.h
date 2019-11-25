
#ifndef FEASST_MONTE_CARLO_ROSENBLUTH_H_
#define FEASST_MONTE_CARLO_ROSENBLUTH_H_

#include "system/include/select_list.h"
#include "system/include/system.h"

namespace feasst {

/**
  Store the energies and Rosenbluth factors necessary for configurational
  bias Monte Carlo trials.
 */
class Rosenbluth {
 public:
  Rosenbluth() {}

  /// Resize the energy, Boltmann factors by the number of steps in the trial.
  void resize(const int num) {
    energy_.resize(num);
    ln_boltzman_.resize(num);
    cumulative_.resize(num);
    stored_.resize(num);
  }

  /// Return the number of steps.
  int num() const { return static_cast<int>(energy_.size()); }

  /// Store the selection for each step.
  void store(const int step, const SelectList& select, const System * system) {
    stored_[step] = select;
    //if (!stored_[step].exchange_indices_positions(select)) {
    //   stored_[step].store(select, system->configuration());
    //}
  }

  /// Set the energy of the step by step.
  void set_energy(const int step, const double energy) {
    energy_[step] = energy;
  }

  /// Compute Boltzmann factors and cumulative probabilities for all steps.
  /// Choose one of the steps based on the probabilities.
  void compute(const double beta, Random * random) {
    for (int step = 0; step < num(); ++step) {
      ln_boltzman_[step] = -beta*energy_[step];
    }
    DEBUG("ln_boltzman " << feasst_str(ln_boltzman_));
    // tot = sum(e^-betaU)
    // ln_tot = ln(sum(e^-betaU)
    // shift by constant, C = -max + 10 to avoid overflow.
    const double shift = 10. - maximum(ln_boltzman_);
    DEBUG("shift " << shift);
    ln_total_rosenbluth_ = 0.;
    for (const double ln : ln_boltzman_) {
      ln_total_rosenbluth_ += exp(ln + shift);
    }
    ln_total_rosenbluth_ = log(ln_total_rosenbluth_) - shift;
    DEBUG("ln_rosen " << ln_total_rosenbluth_);
    DEBUG("energy " << feasst_str(energy_));
    if (ln_total_rosenbluth_ <= -NEAR_INFINITY) {
      chosen_step_ = -1;
      return;
    }
    double accumulator = 0.;
    for (int step = 0; step < num(); ++step) {
      accumulator += exp(ln_boltzman_[step] - ln_total_rosenbluth_);
      cumulative_[step] = accumulator;
    }
    TRACE("cumulative " << feasst_str(cumulative_));
    ASSERT(std::abs(cumulative_.back() - 1.) < 10.*NEAR_ZERO,
      "cumulative probability must end in 1. " <<
      MAX_PRECISION << cumulative_.back());
    TRACE("cumulative " << feasst_str(cumulative_));
    chosen_step_ = random->index_from_cumulative_probability(cumulative_);
  }

  /// Return the stored selection of the step.
  const SelectList& stored(const int step) const { return stored_[step]; }

  /// Return the chosen selection.
  const SelectList& chosen() const {
    ASSERT(chosen_step_ != -1, "error");
    return stored_[chosen_step_];
  }

  /// Return the chosen energy.
  double chosen_energy() const {
    DEBUG("chosen step " << chosen_step_);
    if (chosen_step_ == -1) {
      return energy_[0];
    }
    return energy_[chosen_step_];
  }

  /// Return the natural logrithm of the total rosenbluth factor.
  double ln_total_rosenbluth() const { return ln_total_rosenbluth_; }

  /// Return the energy for a given step.
  double energy(const int step) const { return energy_[step]; }

  /// Return the chosen step
  int chosen_step() const { return chosen_step_; }

  /// Serialize.
  void serialize(std::ostream& ostr) const {
    feasst_serialize_version(507, ostr);
    feasst_serialize(energy_, ostr);
    feasst_serialize(ln_boltzman_, ostr);
    feasst_serialize(cumulative_, ostr);
    feasst_serialize_fstobj(stored_, ostr);
  }

  /// Deserialize
  Rosenbluth(std::istream& istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(version == 507, "version: " << version);
    feasst_deserialize(&energy_, istr);
    feasst_deserialize(&ln_boltzman_, istr);
    feasst_deserialize(&cumulative_, istr);
    feasst_deserialize_fstobj(&stored_, istr);
  }

 private:
  std::vector<double> energy_;
  std::vector<double> ln_boltzman_;
  std::vector<double> cumulative_;
  std::vector<SelectList> stored_;

  // temporary
  double ln_total_rosenbluth_;
  int chosen_step_ = -1;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_ROSENBLUTH_H_

#include <cmath>
#include "utils/include/debug.h"
#include "monte_carlo/include/acceptance.h"

namespace feasst {

double Acceptance::ln_metropolis_prob() const {
  ASSERT(!std::isinf(ln_metropolis_prob_), "ln_metropolis_prob_ is inf");
  return ln_metropolis_prob_;
}

void Acceptance::reset() {
  set_ln_metropolis_prob();
  set_reject();
  set_endpoint();
  energy_new_ = 0.;
  energy_old_ = 0.;
  std::fill(energy_profile_new_.begin(), energy_profile_new_.end(), 0.);
  std::fill(energy_profile_old_.begin(), energy_profile_old_.end(), 0.);
  macrostate_shift_ = 0;
  macrostate_shift_type_ = 0;
  perturbed_.clear();
}

void Acceptance::add_to_perturbed(const Select& select) {
  perturbed_.add(select);
}

void Acceptance::set_perturbed_state(const int state) {
  DEBUG("state " << state);
  perturbed_.set_trial_state(state);
}

void Acceptance::add_to_energy_profile_new(const std::vector<double>& energy) {
  const int current_size = static_cast<int>(energy_profile_new_.size());
  const int new_size = static_cast<int>(energy.size());
  if (current_size < new_size) {
    energy_profile_new_.resize(new_size);
//  } else {
//    ASSERT(current_size == new_size, "err");
  }
  for (int i = 0; i < new_size; ++i) {
    energy_profile_new_[i] += energy[i];
  }
}

void Acceptance::subtract_from_energy_profile_new(const std::vector<double>& energy) {
  const int current_size = static_cast<int>(energy_profile_new_.size());
  const int new_size = static_cast<int>(energy.size());
  if (current_size < new_size) {
    energy_profile_new_.resize(new_size);
//  } else {
//    ASSERT(current_size == new_size, "err");
  }
  for (int i = 0; i < new_size; ++i) {
    energy_profile_new_[i] -= energy[i];
  }
}

void Acceptance::add_to_energy_profile_old(const std::vector<double>& energy) {
  const int current_size = static_cast<int>(energy_profile_old_.size());
  const int new_size = static_cast<int>(energy.size());
  if (current_size < new_size) {
    energy_profile_old_.resize(new_size);
  } else {
    ASSERT(current_size == new_size, "err");
  }
  for (int i = 0; i < new_size; ++i) {
    energy_profile_old_[i] += energy[i];
  }
}

}  // namespace feasst

#include <cmath>
#include "utils/include/debug.h"
#include "utils/include/utils.h"
#include "configuration/include/select.h"
#include "monte_carlo/include/acceptance.h"

namespace feasst {

Acceptance::~Acceptance() {}

double Acceptance::ln_metropolis_prob() const {
  ASSERT(!std::isinf(ln_metropolis_prob_), "ln_metropolis_prob_ is inf");
  return ln_metropolis_prob_;
}

void Acceptance::reset() {
  set_ln_metropolis_prob();
  set_reject();
  set_endpoint();
  energy_new_.resize(1);
  energy_new_[0] = 0.;
  energy_old_.resize(1);
  energy_old_[0] = 0.;
  energy_ref_.resize(1);
  energy_ref_[0] = 0.;
  updated_.resize(1);
  updated_[0] = 0;
  resize(1, 0, &energy_profile_new_);
  fill(0., &energy_profile_new_);
  resize(1, 0, &energy_profile_old_);
  fill(0., &energy_profile_old_);
  macrostate_shift_.resize(1);
  macrostate_shift_[0] = 0;
  macrostate_shift_type_.resize(1);
  macrostate_shift_type_[0] = 0.;
  perturbed_.clear();
  perturbed_.resize(2);  // maximum number of configs
  perturbed_[0] = std::make_shared<Select>();
  perturbed_[1] = std::make_shared<Select>();
}

void Acceptance::add_to_perturbed(const Select& select, const int config) {
  resize_(config, &perturbed_);
  perturbed_[config]->add(select);
}

void Acceptance::set_perturbed_state(const int state, const int config) {
  resize_(config, &perturbed_);
  DEBUG("state " << state);
  perturbed_[config]->set_trial_state(state);
}

const Select& Acceptance::perturbed(const int config) const {
  ASSERT(config < static_cast<int>(perturbed_.size()),
    "config: " << config << " >= size:" << perturbed_.size() <<
    "Consider increasing max num config in reset()");
  return *perturbed_[config];
}

const std::vector<double>& Acceptance::energy_profile_new(
    const int config) const {
  ASSERT(config < static_cast<int>(energy_profile_new_.size()),
    "config:" << config << " >= size:" << energy_profile_new_.size());
  return energy_profile_new_[config];
}

void Acceptance::set_energy_profile_new(const std::vector<double>& energy,
    const int config) {
  if (config == static_cast<int>(energy_profile_new_.size())) {
    energy_profile_new_.resize(energy_profile_new_.size()+1);
  }
  resize_(config, &energy_profile_new_);
  energy_profile_new_[config] = energy;
}

void Acceptance::add_to_energy_profile_new(const std::vector<double>& energy,
    const int config) {
  DEBUG("config " << config);
  resize_(config, &energy_profile_new_);
  const int current_size = static_cast<int>(energy_profile_new_[config].size());
  DEBUG("current_size " << current_size);
  const int new_size = static_cast<int>(energy.size());
  DEBUG("new_size " << new_size);
  if (current_size < new_size) {
    energy_profile_new_[config].resize(new_size);
//  } else {
//    ASSERT(current_size == new_size, "err");
  }
  for (int i = 0; i < new_size; ++i) {
    energy_profile_new_[config][i] += energy[i];
  }
}

void Acceptance::subtract_from_energy_profile_new(
    const std::vector<double>& energy, const int config) {
  resize_(config, &energy_profile_new_);
  DEBUG("config " << config);
  const int current_size = static_cast<int>(energy_profile_new_[config].size());
  DEBUG("current_size " << current_size);
  const int new_size = static_cast<int>(energy.size());
  DEBUG("new_size " << new_size);
  if (current_size < new_size) {
    energy_profile_new_[config].resize(new_size);
//  } else {
//    ASSERT(current_size == new_size, "err");
  }
  for (int i = 0; i < new_size; ++i) {
    energy_profile_new_[config][i] -= energy[i];
  }
}

const std::vector<double>& Acceptance::energy_profile_old(
    const int config) const {
  ASSERT(config < static_cast<int>(energy_profile_old_.size()),
    "config:" << config << " >= size:" << energy_profile_old_.size());
  return energy_profile_old_[config];
}

void Acceptance::set_energy_profile_old(const std::vector<double>& energy,
    const int config) {
  resize_(config, &energy_profile_old_);
  energy_profile_old_[config] = energy;
}

void Acceptance::add_to_energy_profile_old(const std::vector<double>& energy,
    const int config) {
  resize_(config, &energy_profile_old_);
  const int current_size = static_cast<int>(energy_profile_old_[config].size());
  const int old_size = static_cast<int>(energy.size());
  if (current_size < old_size) {
    energy_profile_old_[config].resize(old_size);
//  } else {
//    ASSERT(current_size == old_size, "err");
  }
  for (int i = 0; i < old_size; ++i) {
    energy_profile_old_[config][i] += energy[i];
  }
}

void Acceptance::set_reject(const bool reject) {
  reject_ = reject;
}

double Acceptance::energy_new(const int config) const {
  ASSERT(config < static_cast<int>(energy_new_.size()),
    "config:" << config << " >= size:" << energy_new_.size());
  return energy_new_[config];
}

double Acceptance::energy_ref(const int config) const {
  ASSERT(config < static_cast<int>(energy_ref_.size()),
    "config:" << config << " >= size:" << energy_ref_.size());
  return energy_ref_[config];
}

void Acceptance::set_energy_new(const double energy, const int config) {
  resize_(config, &energy_new_);
  energy_new_[config] = energy;
  resize_(config, &updated_);
  updated_[config] = 1;
}

void Acceptance::set_energy_ref(const double energy, const int config) {
  resize_(config, &energy_ref_);
  energy_ref_[config] = energy;
}

void Acceptance::add_to_energy_new(const double energy, const int config) {
  resize_(config, &energy_new_);
  energy_new_[config] += energy;
  resize_(config, &updated_);
  updated_[config] = 1;
}

double Acceptance::energy_old(const int config) const {
  ASSERT(config < static_cast<int>(energy_old_.size()),
    "config:" << config << " >= size:" << energy_old_.size());
  return energy_old_[config];
}

void Acceptance::set_energy_old(const double energy, const int config) {
  resize_(config, &energy_old_);
  energy_old_[config] = energy;
  resize_(config, &updated_);
  updated_[config] = 1;
}

void Acceptance::add_to_energy_old(const double energy, const int config) {
  resize_(config, &energy_old_);
  energy_old_[config] += energy;
  resize_(config, &updated_);
  updated_[config] = 1;
}

int Acceptance::updated(const int config) const {
  ASSERT(config < static_cast<int>(updated_.size()),
    "config:" << config << " >= size:" << updated_.size());
  return updated_[config];
}

int Acceptance::num_configurations() const {
//  INFO("energy_new_ " << energy_new_.size());
//  INFO("energy_old_ " << energy_old_.size());
//  INFO("updated_ " << updated_.size());
  // return static_cast<int>(energy_new_.size());
  return static_cast<int>(updated_.size());
}

void Acceptance::add_to_macrostate_shift(const int shift, const int config) {
  resize_(config, &macrostate_shift_);
  macrostate_shift_[config] += shift;
}

void Acceptance::set_macrostate_shift_type(const int type, const int config) {
  resize_(config, &macrostate_shift_type_);
  macrostate_shift_type_[config] += type;
}

void Acceptance::sort_perturbed() {
  for (auto pert : perturbed_) {
    pert->sort();
  }
}

}  // namespace feasst

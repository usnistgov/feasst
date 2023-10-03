#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "monte_carlo/include/rosenbluth.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"

namespace feasst {

void Rosenbluth::resize(const int num) {
  energy_.resize(num);
  energy_profile_.resize(num);
  excluded_.resize(num);
  weight_.resize(num);
  cumulative_.resize(num);
  stored_.resize(num);
}

void Rosenbluth::compute(const double beta, Random * random, const bool old) {
  if (num() == 1) {
    chosen_step_ = 0;
    weight_[chosen_step_] =
      -beta*(energy_[chosen_step_] - excluded_[chosen_step_]);
    DEBUG("beta " << beta);
    DEBUG("en " << energy_[chosen_step_]);
    DEBUG("ex " << excluded_[chosen_step_]);
    ln_total_rosenbluth_ = weight_[chosen_step_];
    DEBUG("ln total ros " << ln_total_rosenbluth_);
    return;
  }
  // const double lnk = std::log(num());
  for (int step = 0; step < num(); ++step) {
    weight_[step] = -beta*(energy_[step] - excluded_[step]);
    // weight_[step] = -beta*energy_[step] - lnk;
  }
  DEBUG("num " << num());
  DEBUG("weight " << feasst_str(weight_));
  DEBUG("energy " << feasst_str(energy_));
  DEBUG("excluded " << feasst_str(excluded_));
  // tot = sum(e^-betaU)
  // ln_tot = ln(sum(e^-betaU)
  // shift by constant, C = -max + 10 to avoid overflow.
  const double shift = 10. - maximum(weight_);
  DEBUG("shift " << shift);
  ln_total_rosenbluth_ = 0.;
  for (const double ln : weight_) {
    ln_total_rosenbluth_ += exp(ln + shift);
  }
  ln_total_rosenbluth_ = std::log(ln_total_rosenbluth_) - shift;
  DEBUG("ln_rosen " << ln_total_rosenbluth_);
  DEBUG("energy " << feasst_str(energy_));
  // automatically reject very unlikely acceptances.
  if (ln_total_rosenbluth_ <= -NEAR_INFINITY/pow(10, num())) {
    chosen_step_ = -1;
    return;
  }
  double accumulator = 0.;
  for (int step = 0; step < num(); ++step) {
    accumulator += exp(weight_[step] - ln_total_rosenbluth_);
    cumulative_[step] = accumulator;
  }
  TRACE("cumulative " << feasst_str(cumulative_));
  const double last = cumulative_.back();
  ASSERT(std::abs(last - 1.) < 100000000.*NEAR_ZERO,
    "cumulative probability ends in " << MAX_PRECISION << last <<
    " when it should be 1.");
  for (double& element : cumulative_) {
    element /= last;
  }
  TRACE("cumulative " << feasst_str(cumulative_));
  if (old) {
    chosen_step_ = 0;
  } else {
    chosen_step_ = random->index_from_cumulative_probability(cumulative_);
  }
  TRACE("chosen_step_ " << chosen_step_);
  ln_total_rosenbluth_ -= std::log(num());
}

const Select& Rosenbluth::chosen() const {
  ASSERT(chosen_step_ != -1, "error");
  return stored_[chosen_step_];
}

double Rosenbluth::chosen_energy() const {
  DEBUG("chosen step " << chosen_step_);
  if (chosen_step_ == -1) {
    return energy_[0];
  }
  return energy_[chosen_step_];
}

const std::vector<double>& Rosenbluth::chosen_energy_profile() const {
  DEBUG("chosen step " << chosen_step_);
  if (chosen_step_ == -1) {
    return energy_profile_[0];
  }
  return energy_profile_[chosen_step_];
}

void Rosenbluth::serialize(std::ostream& ostr) const {
  feasst_serialize_version(507, ostr);
  feasst_serialize(energy_, ostr);
  feasst_serialize(energy_profile_, ostr);
  feasst_serialize(excluded_, ostr);
  feasst_serialize(weight_, ostr);
  feasst_serialize(cumulative_, ostr);
  feasst_serialize_fstobj(stored_, ostr);
}

Rosenbluth::Rosenbluth(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 507, "version: " << version);
  feasst_deserialize(&energy_, istr);
  feasst_deserialize(&energy_profile_, istr);
  feasst_deserialize(&excluded_, istr);
  feasst_deserialize(&weight_, istr);
  feasst_deserialize(&cumulative_, istr);
  feasst_deserialize_fstobj(&stored_, istr);
}

void Rosenbluth::set_energy(const int step, const double energy,
    const double excluded) {
  DEBUG("energy: " << energy);
  DEBUG("excluded: " << excluded);
  ASSERT(!std::isinf(energy), "energy: " << energy << " is inf.");
  ASSERT(!std::isnan(energy), "energy: " << energy << " is nan.");
  ASSERT(!std::isinf(excluded), "excluded: " << excluded << " is inf.");
  ASSERT(!std::isnan(excluded), "excluded: " << excluded << " is nan.");
  energy_[step] = energy;
  excluded_[step] = excluded;
}

void Rosenbluth::set_energy_profile(const int step,
  const std::vector<double>& energy) {
  energy_profile_[step] = energy;
}

}  // namespace feasst

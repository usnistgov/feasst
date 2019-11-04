
#include <math.h>
#include <cmath>
#include "math/include/utils_math.h"
#include "utils/include/debug.h"
#include "utils/include/utils_io.h"
#include "flat_histogram/include/ln_probability.h"

namespace feasst {

void LnProbability::normalize() {
  // to avoid overflow in the exp function, shift by the highest value first.
  const double shift = maximum(values_);
  for (double& lnp : values_) {
    lnp -= shift;
  }

  // then normalize by subtracting ln of the sum of the probability
  const double ln_sum = log(sum_probability());
  for (double& lnp : values_) {
    lnp -= ln_sum;
  }
  ASSERT(std::abs(sum_probability() - 1.) < size()*NEAR_ZERO,
    "normalization failure, sum:" << MAX_PRECISION << sum_probability());
}

double LnProbability::sum_probability(const int min,
    const int max) const {
  double sum = 0.;
  ASSERT(min >= 0 && min < size(), "min: " << min << " must be < size: " << size());
  ASSERT(max > min && max < size(), "max: " << max << " must be < size: " << size());
  for (int index = min; index < max + 1; ++index) {
    sum += exp(value(index));
  }
  return sum;
}

void LnProbability::serialize(std::ostream& ostr) const {
  feasst_serialize_version(885, ostr);
  feasst_serialize(values_, ostr);
}

LnProbability::LnProbability(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(885 == version, "mismatch version: " << version);
  feasst_deserialize(&values_, istr);
}

std::vector<int> LnProbability::minima(
    const int num_smooth) const {
  return local_minimum_indices(values(), num_smooth);
}

double LnProbability::saturation_objective(
    const double delta_conjugate,
    int num_smooth) const {
  const std::vector<int> mins = minima();
  const int size = static_cast<int>(mins.size());
  if (size == 0) {
    return pow(values_.front() - values_.back(), 2);
  } else if (size == 1) {
    return saturation_objective_boundary(delta_conjugate, mins[0]);
  } else {
    ERROR("Found " << size << " minimums: " << feasst_str(mins));
  }
}

double LnProbability::saturation_objective_boundary(
    const double delta_conjugate,
    int phase_boundary) const {
  const double prob_low_dens = sum_probability(0, phase_boundary - 1);
  const double prob_high_dens = sum_probability(phase_boundary, size() - 1);
  return pow(log(prob_low_dens) - log(prob_high_dens), 2);
}

}  // namespace feasst

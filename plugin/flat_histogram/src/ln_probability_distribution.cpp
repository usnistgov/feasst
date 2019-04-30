
#include <math.h>
#include <cmath>
#include "core/include/utils_math.h"
#include "core/include/debug.h"
#include "core/include/utils_io.h"
#include "flat_histogram/include/ln_probability_distribution.h"

namespace feasst {

void LnProbabilityDistribution::normalize() {
  // to avoid overflow in the exp function, shift by the highest value first.
  const double shift = maximum(values_);
  for (double& lnp : values_) {
    lnp -= shift;
  }

  // then normalize by subtracting ln of the sum of the probability
  const double ln_sum = log(sum_probability_());
  for (double& lnp : values_) {
    lnp -= ln_sum;
  }
  ASSERT(std::abs(sum_probability_() - 1.) < size()*NEAR_ZERO,
    "normalization failure, sum:" << MAX_PRECISION << sum_probability_());
}

double LnProbabilityDistribution::sum_probability_() {
  double sum = 0.;
  for (double ln_prob : values_) {
    sum += exp(ln_prob);
  }
  return sum;
}

void LnProbabilityDistribution::serialize(std::ostream& ostr) const {
  feasst_serialize_version(885, ostr);
  feasst_serialize(values_, ostr);
}

LnProbabilityDistribution::LnProbabilityDistribution(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(885 == version, "mismatch version: " << version);
  feasst_deserialize(&values_, istr);
}

}  // namespace feasst

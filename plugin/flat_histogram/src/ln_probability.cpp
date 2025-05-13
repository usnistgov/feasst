#include <cmath>
#include <string>
//#include <sstream>
#include <iostream>
#include "utils/include/utils.h"  // is_equal
#include "utils/include/file.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "utils/include/debug.h"
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
  ASSERT(max >= min && max < size(), "max: " << max << " must be < size: "
    << size() << " and >= min " << min);
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

double LnProbability::equilibrium_objective(int num_smooth) const {
  const std::vector<int> mins = minima(num_smooth);
  const int size = static_cast<int>(mins.size());
  if (size == 0) {
    return std::pow(values_.front() - values_.back(), 2);
  } else if (size == 1) {
    return equilibrium_objective_boundary(mins[0]);
  } else {
    FATAL("Found " << size << " minimums: " << feasst_str(mins));
  }
}

double LnProbability::equilibrium_objective_boundary(int phase_boundary) const {
  const double prob_low_dens = sum_probability(0, phase_boundary - 1);
  const double prob_high_dens = sum_probability(phase_boundary, size() - 1);
  return std::pow(log(prob_low_dens) - log(prob_high_dens), 2);
}

bool LnProbability::is_equal(const LnProbability& ln_prob,
    const double tolerance) const {
  if (!feasst::is_equal(values_, ln_prob.values_, tolerance)) {
    INFO("not equal: " << feasst_str(values_));
    INFO("not equal: " << feasst_str(ln_prob.values_));
    return false;
  }
  return true;
}

LnProbability LnProbability::reduce(const int keep_every,
                                    const int shift) const {
  std::vector<double> vals;
  ASSERT(keep_every >= 1, "invalid keep_every: " << keep_every);
  int starting_macro = shift;
  int nattempt = 0;
  while (starting_macro < 0) {
    starting_macro += keep_every;
    ++nattempt;
    ASSERT(nattempt < 1e4, "infinite loop");
  }
  for (int macro = starting_macro; macro < size(); macro += keep_every) {
    vals.push_back(value(macro));
  }
  LnProbability reduced(vals);
  reduced.normalize();
  return reduced;
}

double LnProbability::delta(const int bin) const {
  if (bin == 0) {
    return std::nan("");
  }
  return value(bin) - value(bin - 1);
}

std::vector<double> LnProbability::delta_values() const {
  std::vector<double> dat(size());
  for (int i = 0; i < size(); ++i) {
    dat[i] = delta(i);
  }
  return dat;
}

void LnProbability::read(const std::string& file_name) {
  ASSERT(!file_name.empty(), "no file_name provided");
  ASSERT(file_exists(file_name),
    "file_name:" << file_name << " could not be found.");
  auto file = std::make_shared<std::ifstream>();
  file->open(file_name);
  ASSERT(file->good(), "err");
  std::string line;
  std::getline(*file, line);
  const int nMax = 1e9;
  int i = 0;
  values_.clear();
  values_.push_back(str_to_double(line));
  while ((i != nMax) && (!file->fail())) {
    DEBUG("line: " << line);
    std::getline(*file, line);
    if (!line.empty()) {
      values_.push_back(str_to_double(line));
    }
    ++i;
  }
  normalize();
}

}  // namespace feasst

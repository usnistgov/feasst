#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "flat_histogram/include/window_exponential.h"

namespace feasst {

WindowExponential::WindowExponential(argtype args) : Window(&args) {
  alpha_ = dble("alpha", &args, 1.5);

  DEBUG("Reading optional mins");
  std::string start = "min";
  int index = 0;
  std::stringstream key;
  key << start << index;
  while (used(key.str(), args)) {
    partial_segment_.push_back(dble(key.str(), &args));
    ++index;
    key.str("");
    key << start << index;
    ASSERT(index < 1e8, "index(" << index << ") is very high. Infinite loop?");
  }
  feasst_check_all_used(args);
}

std::vector<double> WindowExponential::segment() const {
  std::vector<double> segment(num() + 1);
  const int num_partial = static_cast<int>(partial_segment_.size());
  DEBUG("num_partial " << num_partial);
  if (num_partial == 0) {
    segment[0] = minimum();
    segment[num()] = maximum();
    const long double diff =
      (std::pow(maximum(), alpha_) - std::pow(segment[0], alpha_))
      /static_cast<double>(num());
    for (int index = 1; index < num(); ++index) {
      segment[index] = std::pow(std::pow(segment[index - 1], alpha_) + diff,
        1./alpha_);
    }
  } else {
    for (int i = 0; i < num_partial; ++i) {
      segment[i] = partial_segment_[i];
    }
    const double last_partial_seg = segment[num_partial-1];
    segment[num()] = maximum();
    const long double diff =
      (std::pow(maximum(), alpha_) - std::pow(last_partial_seg, alpha_))
      /static_cast<double>(num()-num_partial + 1);
    for (int index = num_partial; index < num(); ++index) {
      segment[index] = std::pow(std::pow(segment[index - 1], alpha_) + diff,
        1./alpha_);
    }
  }
  return segment;
}

}  // namespace feasst

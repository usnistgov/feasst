#include <vector>
#include <cmath>
#include "flat_histogram/include/window_exponential.h"

namespace feasst {

WindowExponential::WindowExponential(argtype args) : Window(&args) {
  alpha_ = dble("alpha", &args, 1.5);
  check_all_used(args);
}

std::vector<double> WindowExponential::segment() const {
  std::vector<double> segment(num() + 1);
  segment[0] = minimum();
  segment[num()] = maximum();
  const long double diff =
    (std::pow(maximum(), alpha_) - std::pow(minimum(), alpha_))
    /static_cast<double>(num());
  for (int index = 1; index < num(); ++index) {
    segment[index] = std::pow(std::pow(segment[index - 1], alpha_) + diff,
      1./alpha_);
  }
  return segment;
}

}  // namespace feasst

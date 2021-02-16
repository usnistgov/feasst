#include <vector>
#include <cmath>
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "flat_histogram/include/window.h"
#include "threads/include/thread_omp.h"

namespace feasst {

Window::Window(argtype args) : Window(&args) {
  check_all_used(args);
}
Window::Window(argtype * args) {
  minimum_ = integer("minimum", args, 0);
  maximum_ = integer("maximum", args, -1);
  num_ = integer("num", args, -1);
  if (boolean("num_from_omp", args, false)) {
    ASSERT(num_ == -1, "cannot use both num and num_from_omp args");
    num_ = ThreadOMP().num();
  }
  extra_overlap_ = integer("extra_overlap", args, 0);
  ASSERT(extra_overlap_ >= 0, "extra_overlap: " << extra_overlap_
    << " must be >= 0");
}

//int Window::num() const {
//  ASSERT(num_ > 0, "num: " << num_ << " must be > 0");
//  return num_;
//}

std::vector<std::vector<int> > Window::boundaries() const {
  ASSERT(maximum() > minimum(),
    "maximum " << maximum() << " > minimum: " << minimum());
  std::vector<double> seg = segment();
  std::vector<std::vector<int> > windows(num(), std::vector<int>(2, 0.));
  for (int index = 0; index < num(); ++index) {
    if (index == 0) {
      windows[index][0] = minimum();
    } else {
      windows[index][0] = round(seg[index] - extra_overlap());
    }
    windows[index][1] = round(seg[index + 1]);
  }
  return windows;
}

}  // namespace feasst

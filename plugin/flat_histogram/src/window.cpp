#include <vector>
#include <cmath>
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "flat_histogram/include/window.h"
#include "threads/include/thread_omp.h"

namespace feasst {

Window::Window(const argtype& args) {
  args_.init(args);
  minimum_ = args_.key("minimum").dflt("0").integer();
  maximum_ = args_.key("maximum").dflt("-1").integer();
  num_ = args_.key("num").dflt("-1").integer();
  if (args_.key("num_from_omp").dflt("false").boolean()) {
    ASSERT(num_ == -1, "cannot use both num and num_from_omp args");
    num_ = ThreadOMP().num();
  }
  extra_overlap_ = args_.key("extra_overlap").dflt("0").integer();
  ASSERT(extra_overlap_ >= 0, "extra_overlap: " << extra_overlap_
    << " must be >= 0");
}

int Window::num() const {
  ASSERT(num_ > 0, "num: " << num_ << " must be > 0");
  return num_;
}

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

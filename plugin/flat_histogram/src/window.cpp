#include <vector>
#include <cmath>
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "flat_histogram/include/window.h"
#include "threads/include/thread_omp.h"

namespace feasst {

Window::Window(argtype args) : Window(&args) {
  FEASST_CHECK_ALL_USED(args);
}
Window::Window(argtype * args) {
  minimum_ = integer("minimum", args, 0);
  maximum_ = integer("maximum", args, -1);
  num_ = integer("num", args, -1);
  if (boolean("num_from_omp", args, false)) {
    ASSERT(num_ == -1, "cannot use both num and num_from_omp args");
    num_ = ThreadOMP().num();
  }
  overlap_ = integer("overlap", args, 1);
  ASSERT(overlap_ >= 0, "overlap: " << overlap_ << " must be >= 0");
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
      windows[index][0] = round(seg[index] - overlap() + 1);
    }
    windows[index][1] = round(seg[index + 1]);
  }

  // check that windows are large enough.
  if (overlap() > 0) {
    for (const std::vector<int> win : windows) {
      int max_overlap = overlap() - 1;
      if (win[0] != minimum()) {
        max_overlap *= 2;
      }
      ASSERT(win[1] - win[0] > max_overlap,
        "Windows are assumed to only overlap with the nearest neighbor, " <<
        "but the chosen settings have windows from two or more away " <<
        "overlapping as well. " <<
        "Try changing the input arguments, such as decreasing num, " <<
        "increasing maximum, decreasing minimum, or decreasing overlap. " <<
        "First window of issue: " << feasst_str(win));
    }
  }
  return windows;
}

}  // namespace feasst

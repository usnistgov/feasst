#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "math/include/utils_math.h"
#include "flat_histogram/include/window_custom.h"

namespace feasst {

WindowCustom::WindowCustom(argtype args) : Window(&args) {
  std::string start = "min";
  int index = 0;
  std::vector<double> segment;
  std::stringstream key;
  key << start << index;
  while (used(key.str(), args)) {
    segment.push_back(dble(key.str(), &args));
    ++index;
    key.str("");
    key << start << index;
    ASSERT(index < 1e8, "index(" << index << ") is very high. Infinite loop?");
  }
  segment.push_back(dble("max", &args));
  init_segment_(segment);
  feasst_check_all_used(args);
}

WindowCustom::WindowCustom(const std::vector<double> segment,
                           argtype args) : Window(&args) {
  init_segment_(segment);
  feasst_check_all_used(args);
}

void WindowCustom::init_segment_(std::vector<double> segment) {
  DEBUG("seg " << feasst_str(segment));
  segment_ = segment;
  const int seg_size = static_cast<int>(segment_.size());
  DEBUG("seg size " << seg_size);
  ASSERT(seg_size > 1, "bad size.");
  boundaries_.clear();
  DEBUG("b size " << boundaries_.size());
  for (int win = 0; win < seg_size - 1; ++win) {
    DEBUG("win " << win);
    int extra = -1 + overlap();
    if (win == seg_size - 2) {
      extra = 0;
    }
    boundaries_.push_back({feasst::round(segment_[win]),
                           feasst::round(segment_[win + 1]) + extra});
  }

  // check arguments
  ASSERT(Window::minimum() == 0, "do not use minimum: " << Window::minimum());
  ASSERT(Window::maximum() == -1, "do not use maximum: " << Window::maximum());
  ASSERT(Window::num() == -1, "do not use num: " << Window::num());
}

}  // namespace feasst

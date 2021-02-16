#include "utils/include/debug.h"
#include "flat_histogram/include/window_custom.h"

namespace feasst {

WindowCustom::WindowCustom(const std::vector<double> segment,
                           argtype args) : Window(&args) {
  segment_ = segment;

  // check arguments
  ASSERT(Window::minimum() == 0, "do not use minimum: " << Window::minimum());
  ASSERT(Window::maximum() == -1, "do not use maximum: " << Window::maximum());
  ASSERT(Window::num() == -1, "do not use num: " << Window::num());
  check_all_used(args);
}

}  // namespace feasst

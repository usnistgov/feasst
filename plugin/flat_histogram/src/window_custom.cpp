#include "utils/include/debug.h"
#include "flat_histogram/include/window_custom.h"

namespace feasst {

WindowCustom::WindowCustom(const std::vector<double> segment,
                           const argtype& args) : Window(args) {
  segment_ = segment;

  // check arguments
  Arguments args2(args);
  args2.dont_check();
  if (args2.key("maximum").used()) WARN("maximum argument not used");
  if (args2.key("minimum").used()) WARN("minimum argument not used");
  if (args2.key("num").used()) WARN("num argument not used");
}

}  // namespace feasst

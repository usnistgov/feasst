#include "utils/include/utils.h"

namespace feasst {

bool is_equal(double val1, double val2, const double tolerance) {
  if (std::abs(val1 - val2) > tolerance) {
    return false;
  }
  return true;
}

}  // namespace feasst

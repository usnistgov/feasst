#include "utils/include/utils.h"
//#include "utils/include/debug.h"

namespace feasst {

bool is_equal(double val1, double val2, const double tolerance) {
  if (std::abs(val1 - val2) > tolerance) {
    return false;
  }
  return true;
}

bool is_equal_fixed_tolerance(const float val1, const float val2) {
  if (val1 == val2) {
    return true;
  }
  //INFO("val1: " << val1 << " val2: " << val2);
  return false;
}

}  // namespace feasst

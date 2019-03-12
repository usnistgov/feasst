
#include <cmath>
#include "core/include/debug.h"
#include "core/include/tunable.h"

namespace feasst {

Tunable::Tunable() {
  set_percent_change();
  set_target();
}

void Tunable::set_min_and_max(const double min, const double max) {
  is_bound_ = true;
  min_ = min;
  max_ = max;
}
  
void Tunable::set_value(const double value) {
  if (is_bound_) {
    if (value > min_ && value < max_) {
      value_ = value;
    }
  } else {
    value_ = value;
  }
}

void Tunable::set_percent_change(const double percent) {
  ASSERT(std::abs(percent) < 1., "|percent| as decimal should be less than 1");
  percent_change_ = percent;
}
  
void Tunable::tune(const double actual) {
  if (is_enabled_) {
    double value = value_;
    if (actual < target_) {
      value *= 1. - percent_change_;
    } else {
      value *= 1. + percent_change_;
    }
    set_value(value);
  }
}

}  // namespace feasst

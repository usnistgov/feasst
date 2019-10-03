
#include <cmath>
#include "utils/include/debug.h"
#include "monte_carlo/include/tunable.h"
#include "utils/include/utils_io.h"
#include "math/include/constants.h"

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
  DEBUG("setting val: " << value);
  if (is_bound_) {
    if (value > min_ and value < max_) {
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
    DEBUG("target: " << target_ << " actual: " << actual);
    if (actual < target_) {
      value *= 1. - percent_change_;
    } else {
      value *= 1. + percent_change_;
    }
    set_value(value);
  }
}


void Tunable::serialize(std::ostream& ostr) const {
  feasst_serialize_version(265, ostr);
  feasst_serialize(is_enabled_, ostr);
  feasst_serialize(value_, ostr);
  feasst_serialize(is_bound_, ostr);
  feasst_serialize(max_, ostr);
  feasst_serialize(min_, ostr);
  feasst_serialize(target_, ostr);
  feasst_serialize(percent_change_, ostr);
}

Tunable::Tunable(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 265, "version");
  feasst_deserialize(&is_enabled_, istr);
  feasst_deserialize(&value_, istr);
  feasst_deserialize(&is_bound_, istr);
  feasst_deserialize(&max_, istr);
  feasst_deserialize(&min_, istr);
  feasst_deserialize(&target_, istr);
  feasst_deserialize(&percent_change_, istr);
}

bool Tunable::is_equal(const Tunable& tunable) const {
  if (is_enabled_ != tunable.is_enabled_) {
    DEBUG("unequal is_enabled:" << is_enabled_ << " " << tunable.is_enabled_);
    return false;
  }
  if (std::abs(value_ - tunable.value_) > NEAR_ZERO) {
    DEBUG("unequal value:" << value_ << " " << tunable.value_);
    return false;
  }
  return true;
}

}  // namespace feasst

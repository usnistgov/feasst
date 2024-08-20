
#include <cmath>
#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/tunable.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"

namespace feasst {

Tunable::Tunable(argtype * args) {
  set_value(dble("tunable_param", args, 0.1));
  set_target(dble("tunable_target_acceptance", args, 0.25));
  set_percent_change(dble("tunable_percent_change", args, 1));
}
Tunable::Tunable(argtype args) : Tunable(&args) {
  feasst_check_all_used(args);
}

void Tunable::set_min_and_max(const double min, const double max) {
  is_bound_ = true;
  min_ = min;
  max_ = max;
}

void Tunable::set_value(const double value) {
  DEBUG("setting val: " << value);
  if (is_bound_) {
    if (value > min_ && value < max_) {
      value_ = value;
    }
  } else {
    value_ = value;
  }
}

void Tunable::set_percent_change(const double percent) {
  percent_change_ = percent;
}

void Tunable::tune(const double actual) {
  DEBUG("is_enabled " << is_enabled_);
  if (is_enabled_) {
    double value = value_;
    DEBUG("target: " << target_ << " actual: " << actual);
    DEBUG("old value " << value);
    value *= 1 + percent_change_*(actual - target_);
    DEBUG("new value " << value);
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

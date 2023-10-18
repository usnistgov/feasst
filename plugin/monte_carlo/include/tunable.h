
#ifndef FEASST_MONTE_CARLO_TUNABLE_H_
#define FEASST_MONTE_CARLO_TUNABLE_H_

#include "utils/include/arguments.h"

namespace feasst {

/**
  A tunable parameter attempts to reach a target value of acceptance by
  modifying it by a percentage.

  The parameter value is changed by the following formula:

  value *= 1 + percent_change*(actual - target)

  where actual is the average trial acceptance and target is the target trial
  acceptance.
  The relationship between the target and the parameter is assumed by the
  sign of percent_change.
  By default, the precent_change = 1, which leads to an inverse relationship
  between the parameter value and the trial acceptance.
  For example, a maximum displacement parameter with an actual acceptance of 32\%
  and a target acceptance of 25\% would be increased by 7\% in an attempt to
  bring the acceptance down to the target value.
 */
class Tunable {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - tunable_param: initial value of the tunable parameter (default: 0.1).
    - tunable_target_acceptance: optionally set target acceptance (default: 0.25).
    - tunable_percent_change: optionally set the percent change (default: 1).
   */
  explicit Tunable(argtype args = argtype());
  explicit Tunable(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  double target() const { return target_; }
  double percent_change() const { return percent_change_; }
  bool is_bound() const { return is_bound_; }

  /// By default, the value is not bound.
  /// This function sets the minimum and maximum bounds.
  void set_min_and_max(const double min, const double max);

  /// Return the above parameter.
  double max() const { return max_; }

  /// Return the above parameter.
  double min() const { return min_; }

  /// Set the value, optionally subject to bounds.
  void set_value(const double value);

  /// Return the value.
  double value() const { return value_; }

  /// Disable tuning.
  void disable() { is_enabled_ = false; }

  /// Enable tuning.
  void enable() { is_enabled_ = true; }

  /// Return if enabled.
  bool is_enabled() const { return is_enabled_; }

  /// A positive percent change indicates that the value of the parameter and
  /// the target have an inverse relationship.
  /// Thus, if the actual value is greater than the target, the value will be
  /// increased upon tuning.
  void set_percent_change(const double percent = 1);

  /// Set the target to reach for tuning.
  void set_target(const double target = 0.25) { target_ = target; }

  /// Change the value to attempt to reach a point where the actual is equal to
  /// the target.
  void tune(const double actual);

  // Check if approximately equal to given tunable.
  bool is_equal(const Tunable& tunable) const;

  /// Serialize.
  void serialize(std::ostream& ostr) const;

  /// Construct from serialization.
  explicit Tunable(std::istream& istr);

  //@}
 private:
  bool is_enabled_ = true;
  double value_ = 0.;
  bool is_bound_ = false;
  double max_ = 0.;
  double min_ = 0.;
  double target_;
  double percent_change_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TUNABLE_H_

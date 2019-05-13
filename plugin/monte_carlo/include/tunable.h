
#ifndef FEASST_MONTE_CARLO_TUNABLE_H_
#define FEASST_MONTE_CARLO_TUNABLE_H_

namespace feasst {

/**
  A tunable parameter attempts to reach a target value of acceptance (or
  generic measure) by modifying it by a percentage.
  The relationship between the target and the parameter is assumed by the
  sign of the percent change (by default, inverse relationship assumed).
 */
class Tunable {
 public:
  Tunable();

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

  /// A positive percent change indicates that the value of the parameter and
  /// the target have an inverse relationship.
  /// Thus, if the actual value is greater than the target, the value will be
  /// increased upon tuning.
  void set_percent_change(const double percent = 0.05);

  /// Set the target to reach for tuning.
  void set_target(const double target = 0.25) { target_ = target; }

  /// Change the value to attempt to reach a point where the actual is equal to
  /// the target.
  void tune(const double actual);

 private:
  bool is_enabled_ = true;
  double value_ = 0.;
  bool is_bound_ = false;
  double max_;
  double min_;
  double target_;
  double percent_change_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TUNABLE_H_

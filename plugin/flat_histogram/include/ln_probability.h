
#ifndef FEASST_FLAT_HISTOGRAM_LN_PROBABILITY_H_
#define FEASST_FLAT_HISTOGRAM_LN_PROBABILITY_H_

#include <vector>

namespace feasst {

/**
  Natural logarithm of a probability distribution.
 */
class LnProbability {
 public:
  LnProbability() {}

  /// Construct with given values.
  LnProbability(const std::vector<double>& values) {
    values_ = values; }

  /// Return the value of the bin.
  double value(const int bin) const { return values_[bin]; }

  /// Set the value of the bin.
  void set_value(const int bin, const double value) { values_[bin] = value; }

  /// Add to the value of the bin.
  void add(const int bin, const double value_to_add) {
    values_[bin] += value_to_add; }

  /// Set the number of macrostates.
  void resize(const int size) { values_.resize(size); }

  /// Return the values.
  const std::vector<double> values() const { return values_; }

  /// Return the size of the distribution.
  int size() const { return static_cast<int>(values_.size()); }

  /// Return the sum of the probability from min to max indices.
  double sum_probability(const int min, const int max) const;

  /// Return the sum of the probability of all values.
  double sum_probability() const { return sum_probability(0, size() - 1); }

  /// Normalize such that the sum of the probability is unity.
  void normalize();

  /** Return the local minimum indices which are canditates for phase boundary.
      The returned indices must be global minimum between +/- num_smooth. */
  std::vector<int> minima(
    const int num_smooth = 10) const;

  /** Return the objective function to minimize in order to obtain saturation
      conditions.
      If no phase boundaries, return the squared difference between
      the minimum and maximum values.
      If more than one, no current implementation.
      If exactly one, use objective function below. */
  double saturation_objective(const double delta_conjugate,
    /// Number of macrostates which define 'local' region for finding the
    /// minimum.
    const int num_smooth = 10) const;

  /** Return the objective function to minimize in order to obtain saturation
      conditions.
      The objective is the squared difference between the natural logs of the
      probabilities of the two phases.
    */
  double saturation_objective_boundary(
    /// See FlatHistogram::reweight()
    const double delta_conjugate,
    /// Index which is the boundary between phases.
    /// The exact index value is included in the second phase.
    int phase_boundary) const;

  bool is_equal(const LnProbability& ln_prob, const double tolerance) const;

  void serialize(std::ostream& ostr) const;
  LnProbability(std::istream& istr);

 private:
  std::vector<double> values_;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_LN_PROBABILITY_H_


#ifndef FEASST_FLAT_HISTOGRAM_LN_PROBABILITY_DISTRIBUTION_H_
#define FEASST_FLAT_HISTOGRAM_LN_PROBABILITY_DISTRIBUTION_H_

#include <vector>

namespace feasst {

/**
  Natural logarithm of a probability distribution.
 */
class LnProbabilityDistribution {
 public:
  LnProbabilityDistribution() {}

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

  /// Normalize such that the sum of the probability is unity.
  void normalize();

  void serialize(std::ostream& ostr) const;

  LnProbabilityDistribution(std::istream& istr);

 private:
  std::vector<double> values_;

  // compute the sum of the probability
  double sum_probability_();
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_LN_PROBABILITY_DISTRIBUTION_H_


#ifndef FEASST_MATH_HISTOGRAM_H_
#define FEASST_MATH_HISTOGRAM_H_

#include <deque>
#include <memory>
#include "math/include/formula.h"
#include "utils/include/arguments.h"

namespace feasst {

/**
  A one-dimensional histogram.
  The histogram is defined by a contiguous set of "bins" which may be of
  constant or variable width.
  Each time a value is added to the histogram, the bin is found and incremented.
 */
class Histogram {
 public:
  /**
    args:
    - width: constant bin width
    - max: maximum value
    - min: minimum value (default: 0)
   */
  Histogram(const argtype& args = argtype());

  /**
    Depending on the way the bins are defined, the histogram may be told how to
    automatically expand in size if the value is outside of the current bounds.
    Otherwise, adding a value outside of the current bounds will result in an
    error.
   */
  bool expandable() const { return expandable_; }

  /// Set the size of the bins according to a formula, where the initial bin is
  /// defined by the formula evaluated at zero.
  /// This makes the histogram expandable.
  void set_bin_size(const std::shared_ptr<Formula> bin_size);

  /// Set a constant bin width and fix the center of the initial bin.
  /// This makes the histogram expandable.
  void set_width_center(const double width, const double center);

  /// Set the edges or boundaries manually.
  /// This makes the histogram not expandable.
  void set_edges(const std::deque<double> edges);

  /// Force the histogram not to be expandable.
  void set_not_expandable() { expandable_ = false; }

  /// The edges or boundaries of the contiguous set of bins.
  std::deque<double> edges() const { return edges_; }

  /// Return the maximum edge.
  double max() const { return edges_.back(); }

  /// Return the minimum edge.
  double min() const { return edges_.front(); }

  /// Return the number of bins.
  int size() const { return histogram_.size(); }

  /// Return the bin for a given value.
  int bin(const double value) const;

  /// Return the center of the bin (e.g., the average of the boundaries).
  double center_of_bin(const int bin) const;

  /// Return the center of the last bin.
  double center_of_last_bin() const {
    return center_of_bin(size() - 1);
  }

  /// Update the histogram by adding a value.
  void add(const double value);

  /// Return the histogram.
  std::deque<double> histogram() const { return histogram_; }

  void serialize(std::ostream& ostr) const;
  Histogram(std::istream& istr);

 private:
  std::deque<double> histogram_;
  std::deque<double> edges_;
  bool expandable_ = false;
  std::shared_ptr<Formula> bin_size_;

  // optimization for constant width
  int is_constant_width_ = 0;

  Arguments args_;

  void set_expandable_() { expandable_ = true; }
};

}  // namespace feasst

#endif  // FEASST_MATH_HISTOGRAM_H_


#ifndef FEASST_FLAT_HISTOGRAM_WINDOW_H_
#define FEASST_FLAT_HISTOGRAM_WINDOW_H_

#include <vector>
#include "utils/include/arguments.h"

namespace feasst {

/**
  Macrostate ranges are often broken into separate windows for parallelization.
  The relative size of these windows can be optimized depending upon the system.
  Each window must share atleast one macrostate with its neighbors in order to
  recover the free energy of the entire range.

  A segment is a countinuous, theoretical formula for analytically spacing
  macrostate ranges.
  This line segment includes both the minimum, maximum, and intermediate values.
  Thus, a 4 window segment would contain 5 values.

  Boundaries are integer macrostate minimums and maximums utilized by each
  window.
  Thus, 4 window boundaries would contain 4 min and max values for each window.
  Boundaries use the given line segment, either from formula or user input,
  and then account for rounding and extra overlap between windows.

  For example, 4 windows with segments given by [0, 100.00, 141.42, 173.21, 200]
  with 4 extra overlap and rounding would result in the following boundaries:
  [[0, 100], [96, 141], [137, 173], [169, 200]] where each of the 4 window
  is expressed as a [min, max] pair.

  Note that TransitionMatrix should have atleast 3 overlap because the first
  and last macrostates in the window automatically reject deletions/insertions,
  respectively, which could skew the calculated probability distributions.
  WangLandau does not have this issue.
 */
class Window {
 public:
  /**
    args:
    - minimum: minimum in macrostate range (default: 0).
    - maximum: maximum in macrostate range.
    - num: number of windows
    - extra_overlap: additional overlap of the windows in addition to the one
      macrostate that they must share (default: 0).
   */
  explicit Window(const argtype& args = argtype());

  /// Return the minimum.
  int minimum() const { return minimum_; }

  /// Return the maximum.
  int maximum() const { return maximum_; }

  /// Return the number of windows.
  int num() const { return num_; }

  /// Return the extra overlap.
  int extra_overlap() const { return extra_overlap_; }

  /// Return the continuous, segmented boundaries of the range.
  /// This should be return num + 1 boundaries, to include global min and max.
  virtual std::vector<double> segment() const = 0;

  /// Return the boundaries as a vector of vectors.
  std::vector<std::vector<int> > boundaries() const;

 protected:
  Arguments args_;

 private:
  int minimum_, maximum_, num_, extra_overlap_;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WINDOW_H_

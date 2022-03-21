
#ifndef FEASST_FLAT_HISTOGRAM_WINDOW_H_
#define FEASST_FLAT_HISTOGRAM_WINDOW_H_

#include <vector>
#include "utils/include/arguments.h"

namespace feasst {

/**
  Macrostate ranges are often broken into separate windows for parallelization.
  The relative size of these windows can be optimized depending upon the system.

  A segment is a countinuous, theoretical formula for analytically spacing
  macrostate ranges.
  This line segment includes both the minimum, maximum, and intermediate values.
  Thus, a 4 window segment would contain 5 values.

  Boundaries are integer macrostate minimums and maximums utilized by each
  window.
  Thus, 4 window boundaries would contain 4 min and max values.
  Boundaries use the given line segment, either from formula or user input,
  and then account for rounding and extra overlap between windows.

  For example, 4 windows with segments given by [0, 100.00, 141.42, 173.21, 200]
  with 5 overlap and rounding would result in the following boundaries:
  [[0, 100], [96, 141], [137, 173], [169, 200]] where each of the 4 window
  is expressed as a [min, max] pair.

  In order to recover the free energy over all of the windows, it is possible
  to splice either (1) the probability distribution or (2) the collection matrix
  elements.

  For (1) probability distribution splicing, each window must share atleast one
  macrostate (overlap) with its neighbors.

  For (2) collection matrix splicing, windows cannot share any macrostates but
  the boundaries must be adjacent and all trial moves of the windows must be
  exactly the same.
 */
class Window {
 public:
  /**
    args:
    - minimum: minimum in macrostate range (default: 0).
    - maximum: maximum in macrostate range (default: -1).
    - num: number of windows (default: -1).
    - num_from_omp: obtain num from OMP threads (default: false).
    - overlap: number of macrostate overlaps between windows (default: 1).
   */
  explicit Window(argtype args = argtype());
  explicit Window(argtype * args);

  /// Return the minimum.
  virtual int minimum() const { return minimum_; }

  /// Return the maximum.
  virtual int maximum() const { return maximum_; }

  /// Return the number of windows.
  virtual int num() const { return num_; }

  /// Return the overlap.
  int overlap() const { return overlap_; }

  /// Return the continuous, segmented boundaries of the range.
  /// This should be return num + 1 boundaries, to include global min and max.
  virtual std::vector<double> segment() const = 0;

  /// Return the boundaries as a vector of vectors.
  std::vector<std::vector<int> > boundaries() const;

  virtual ~Window() {}

 private:
  int minimum_, maximum_, num_, overlap_;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WINDOW_H_

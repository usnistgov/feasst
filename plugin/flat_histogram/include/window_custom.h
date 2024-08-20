
#ifndef FEASST_FLAT_HISTOGRAM_WINDOW_CUSTOM_H_
#define FEASST_FLAT_HISTOGRAM_WINDOW_CUSTOM_H_

#include <vector>
#include "flat_histogram/include/window.h"

namespace feasst {

/**
  Input custom window values directly by manual input of the segment, defined
  in Window.
 */
class WindowCustom : public Window {
 public:
  //@{
  /** @name Arguments
    - min[i]: minimum macrostate value of the i-th window.
    - max: maximum of the largest window.
    - Window arguments.
   */
  explicit WindowCustom(argtype args = argtype());

  //@}
  /** @name Public Functions
   */
  //@{

  /// Construct by manual input of the segments, defined in Window.
  WindowCustom(const std::vector<double> segment,
               argtype args = argtype());
  std::vector<double> segment() const override { return segment_; }
  std::vector<std::vector<int> > boundaries() const override {
    return boundaries_; }
  int minimum() const override { return static_cast<int>(segment_.front()); }
  int maximum() const override { return static_cast<int>(segment_.back()); }
  int num() const override { return static_cast<int>(segment_.size() - 1); }
  virtual ~WindowCustom() {}

  //@}
 private:
  std::vector<double> segment_;
  std::vector<std::vector<int> > boundaries_;
  void init_segment_(std::vector<double> segment);
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WINDOW_CUSTOM_H_

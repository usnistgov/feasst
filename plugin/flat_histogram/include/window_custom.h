
#ifndef FEASST_FLAT_HISTOGRAM_WINDOW_CUSTOM_H_
#define FEASST_FLAT_HISTOGRAM_WINDOW_CUSTOM_H_

#include <vector>
#include "utils/include/arguments.h"
#include "flat_histogram/include/window.h"

namespace feasst {

/**
  Input custom window values directly by manual input of the segment, defined
  in Window.
 */
class WindowCustom : public Window {
 public:
  /// Construct by manual input of the segments, defined in Window.
  WindowCustom(const std::vector<double> segment,
               const argtype& args = argtype());

  std::vector<double> segment() const override { return segment_; }
  int minimum() const override { return static_cast<int>(segment_.front()); }
  int maximum() const override { return static_cast<int>(segment_.back()); }
  int num() const override { return static_cast<int>(segment_.size() - 1); }
  virtual ~WindowCustom() {}

 private:
  std::vector<double> segment_;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WINDOW_CUSTOM_H_

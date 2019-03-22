
#ifndef FEASST_FLAT_HISTOGRAM_MACROSTATE_H_
#define FEASST_FLAT_HISTOGRAM_MACROSTATE_H_

#include "core/include/system.h"
#include "core/include/criteria.h"
#include "core/include/histogram.h"

namespace feasst {

/**
  The macrostate is the statistical mechanical variable to which we apply bias
  using flat-histogram methods.
  Thus, for example, in the grand canonical ensemble with a single component,
  the macrostate is the number of particles.
  To apply the flat histogram methods, the macrostate is be broken into a
  contiguous series of "bins".
 */
class Macrostate {
 public:
  Macrostate(const Histogram& histogram) {
    set(histogram);
  }

  /// Return the current value of the macrostate.
  virtual double value(const System* system, const Criteria* criteria) = 0;

  /// Set the bins of the macrostate by providing a Histogram.
  /// This is required before the macrostate can be used for flat histogram
  /// methods.
  /// The histogram only serves to determine the bins, and should not be
  /// expanded or have values added during the course of the simulation.
  void set(const Histogram histogram) { histogram_ = histogram; }

  /// Return the histogram.
  const Histogram& histogram() const { return histogram_; }

  /// Return the current bin of the macrostate.
  int bin(const System* system, const Criteria* criteria) {
    return histogram_.bin(value(system, criteria)); }

  /// Return whether the current system macrostate is within permissible range
  /// given by the input histogram.
  bool is_in_range(const System* system, const Criteria* criteria);

  virtual ~Macrostate() {}

 private:
  Histogram histogram_;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_MACROSTATE_H_

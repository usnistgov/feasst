
#ifndef FEASST_FLAT_HISTOGRAM_MACROSTATE_H_
#define FEASST_FLAT_HISTOGRAM_MACROSTATE_H_

#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "math/include/histogram.h"

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

  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Macrostate> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Macrostate> >& deserialize_map();
  std::shared_ptr<Macrostate> deserialize(std::istream& istr);
  Macrostate(std::istream& istr);
  virtual ~Macrostate() {}

 protected:
  void serialize_macrostate_(std::ostream& ostr) const;

 private:
  Histogram histogram_;
};

/// Segment an range into pieces by exponential scaling.
inline std::vector<double> segment(
    const double min,   //!< minimum in range
    const double max,   //!< maximum in range
    const int num,      //!< number of segments
    const double exp    //!< exponential parameter
    ) {
  ASSERT(num > 0, "num(" << num << ") must be > 0.");
  std::vector<double> segment(num + 1);
  segment[0] = min;
  segment[num] = max;
  const long double exp_diff =
    (pow(max, exp) - pow(min, exp))/static_cast<double>(num);
  for (int index = 1; index < num; ++index) {
    segment[index] = pow(pow(segment[index - 1], exp) + exp_diff, 1./exp);
  }
  return segment;
}

/// Segment a range into windows by exponential scaling.
inline std::vector<std::vector<int> > window(
    const int min,   //!< minimum in range
    const int max,   //!< maximum in range
    const int num,      //!< number of segments
    const double exp,    //!< exponential parameter
    const int extra_overlap = 0 //!< extra overlap between windows
    ) {
  std::vector<double> boundaries = segment(min, max, num, exp);
  std::vector<std::vector<int> > windows(num, std::vector<int>(2, 0.));
  for (int index = 0; index < num; ++index) {
    if (index == 0) {
      windows[index][0] = min;
    } else {
      windows[index][0] = round(boundaries[index] - extra_overlap);
    }
    windows[index][1] = round(boundaries[index + 1]);
  }
  return windows;
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_MACROSTATE_H_

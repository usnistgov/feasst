
#ifndef FEASST_FLAT_HISTOGRAM_MACROSTATE_H_
#define FEASST_FLAT_HISTOGRAM_MACROSTATE_H_

#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "math/include/histogram.h"
#include "flat_histogram/include/constraint.h"

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
  Macrostate(const Histogram& histogram,
    /**
      soft_max : optionally, set a soft maximum (default: same as histogram).
                 These soft limits may be changed during a simulation.
      soft_min : minimum as described above (default: same as histogram).
     */
    const argtype& args = argtype()) {
    set(histogram);

    // soft limits
    Arguments args_(args);
    is_soft_bound_ = false;
    soft_min_ = soft_max_ = 0.;
    if (args_.key("soft_max").used()) {
      is_soft_bound_ = true;
      soft_max_ = args_.integer();
      if (args_.key("soft_min").used()) {
        soft_min_ = args_.integer();
      } else {
        soft_min_ = histogram_.min();
      }
    }
  }

  /// Set the bins of the macrostate by providing a Histogram.
  /// This is required before the macrostate can be used for flat histogram
  /// methods.
  /// The histogram only serves to determine the bins, and should not be
  /// expanded or have values added during the course of the simulation.
  void set(const Histogram histogram) { histogram_ = histogram; }

  /// Return the histogram.
  const Histogram& histogram() const { return histogram_; }

  /// Return the soft maximum (default: histogram max).
  const int soft_max() const;

  /// Return the soft minimum (default: histogram min).
  const int soft_min() const;

  /// Return the current value of the macrostate.
  virtual double value(const System* system, const Criteria* criteria) = 0;

  /// Return the current bin of the macrostate.
  int bin(const System* system, const Criteria* criteria) {
    return histogram_.bin(value(system, criteria)); }

  /// Return the value of the bin.
  double value(const int bin) const { return histogram_.center_of_bin(bin); }

  /// Add a constraint.
  void add(std::shared_ptr<Constraint> constraint) {
    constraints_.push_back(constraint); }

  /// Return whether the current system macrostate is within permissible range
  /// given by the input histogram and check any additional constraints.
  bool is_allowed(const System* system, const Criteria* criteria);

  /// Swap the soft bounds with another macrostate.
  void swap_soft_bounds(Macrostate * macrostate);

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
  bool is_soft_bound_;
  int soft_max_;
  int soft_min_;
  std::vector<std::shared_ptr<Constraint> > constraints_;
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

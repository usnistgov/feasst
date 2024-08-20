
#ifndef FEASST_FLAT_HISTOGRAM_MACROSTATE_H_
#define FEASST_FLAT_HISTOGRAM_MACROSTATE_H_

#include <map>
#include <string>
#include <memory>

namespace feasst {

class Acceptance;
class Criteria;
class Histogram;
class System;

typedef std::map<std::string, std::string> argtype;

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
  //@{
  /** @name Arguments
    - soft_macro_max : optionally, set a soft maximum (default: last histogram bin).
      These soft limits may be changed during a simulation.
      Note that this max is a macrostate value, not an integer bin index.
    - soft_macro_min : minimum as described above (default: same as histogram).
      This argument also requires the use of soft_macro_max.
    - Histogram arguments.
   */
  explicit Macrostate(argtype args = argtype());

  //@}
  /** @name Public Functions
   */
  //@{

  /// Arguments as described above, but with explicit histogram object.
  explicit Macrostate(const Histogram& histogram, argtype args = argtype());
  Macrostate(const Histogram& histogram, argtype * args);

  /**
    Set the bins of the macrostate by providing a Histogram.
    This is required before the macrostate can be used for flat histogram
    methods.
    The histogram only serves to determine the bins, and should not be
    expanded or have values added during the course of the simulation.
   */
  void set(const Histogram histogram);

  /// Return the histogram.
  const Histogram& histogram() const;

  /// Return the soft maximum as an integer bin index, not a macrostate.
  const int soft_max() const { return soft_max_; }

  /// Return the soft minimum as an integer bin index, not a macrostate.
  const int soft_min() const { return soft_min_; }

  /// Return the number of macrostates in the soft range.
  const int num_macrostates_in_soft_range() const {
    return soft_max_ - soft_min_ + 1; }

  /// Return the current value of the macrostate.
  virtual double value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const = 0;

  /// Return the current bin of the macrostate.
  int bin(const System& system,
      const Criteria& criteria,
      const Acceptance& acceptance) const;

  /// Return the value of the bin.
  double value(const int bin) const;

  /// Return whether the current system macrostate is within permissible range
  /// given by the input histogram.
  bool is_allowed(const System& system,
                  const Criteria& criteria,
                  const Acceptance& acceptance) const;

//  // Swap the soft bounds with another macrostate.
//  void swap_soft_bounds(Macrostate * macrostate);

  // HWH hackish adjust_bounds interface. See CollectionMatrixSplice.
  int set_soft_max(const int index, const System& sys,
                   const Criteria& criteria);
  int set_soft_min(const int index, const System& sys,
                   const Criteria& criteria);
  void add_to_soft_max(const int num) { soft_max_ += num; }
  void remove_from_soft_min(const int num) { soft_min_ -= num; }

  const std::string& class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Macrostate> create(std::istream& istr) const;
  virtual std::shared_ptr<Macrostate> create(argtype * args) const;
  std::map<std::string, std::shared_ptr<Macrostate> >& deserialize_map();
  std::shared_ptr<Macrostate> deserialize(std::istream& istr);
  std::shared_ptr<Macrostate> factory(const std::string name, argtype * args);
  explicit Macrostate(std::istream& istr);
  virtual ~Macrostate() {}

  //@}
 protected:
  void serialize_macrostate_(std::ostream& ostr) const;
  std::string class_name_;

 private:
  std::unique_ptr<Histogram> histogram_;
  int soft_max_;
  int soft_min_;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_MACROSTATE_H_

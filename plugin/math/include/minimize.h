
#ifndef FEASST_MATH_MINIMIZE_H_
#define FEASST_MATH_MINIMIZE_H_

#include <map>
#include <memory>
#include <string>
#include <sstream>
#include "utils/include/arguments.h"

namespace feasst {

class Formula;

/**
  Minimize the objective function within the bounds up to some tolerance.
 */
class Minimize {
 public:
  /**
    args:
    - tolerance: solve within this tolerance.
    - lower: lower bound.
    - upper: upper bound.
   */
  Minimize(const argtype& args = argtype());

  /// Return the tolerance
  double tolerance() const { return tolerance_; }

  /// Return lower bound
  double lower() const { return lower_; }

  /// Set lower bound
  void set_lower(const double lower) { lower_ = lower; }

  /// Return upper bound
  double upper() const { return upper_; }

  /// Set upper bound
  void set_upper(const double upper) { upper_ = upper; }

  /// Return reduced bounds containing a minimum.
  virtual void bracket(double * lower, double * upper, Formula * formula) = 0;

  /// Return the minimum of the Formula.
  double minimum(Formula * formula);

  // serialization
  virtual std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Minimize> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Minimize> >& deserialize_map();
  std::shared_ptr<Minimize> deserialize(std::istream& istr);
  explicit Minimize(std::istream& istr);
  virtual ~Minimize() {}

 protected:
  std::string class_name_ = "Minimize";
  void serialize_solver_(std::ostream& ostr) const;
  Arguments args_;

 private:
  double tolerance_;
  double lower_ = 0;
  double upper_ = 0;
};

}  // namespace feasst

#endif  // FEASST_MATH_MINIMIZE_H_

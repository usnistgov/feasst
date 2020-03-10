
#ifndef FEASST_MATH_FORMULA_H_
#define FEASST_MATH_FORMULA_H_

#include <map>
#include <memory>
#include <string>
#include <sstream>
#include "utils/include/arguments.h"

namespace feasst {

/**
  Define a formula, used in this code for an expandable histogram.
 */
/// HWH Formula remains to be well tested, also in context of histogram.
class Formula {
 public:
  Formula(
    /**
      x0 : central reference position (default: 0).
     */
    const argtype& args = argtype()) {
    args_.init(args);
    set_x0(args_.key("x0").dflt("0").dble());
  }
  void set_x0(const double x0) { x0_ = x0; }
  double x0() const { return x0_; }
  virtual double evaluate(const double x) = 0;
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Formula> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Formula> >& deserialize_map();
  std::shared_ptr<Formula> deserialize(std::istream& istr);
  explicit Formula(std::istream& istr);
  virtual ~Formula() {}

 protected:
  void serialize_formula_(std::ostream& ostr) const;
  Arguments args_;

 private:
  double x0_;
};

}  // namespace feasst

#endif  // FEASST_MATH_FORMULA_H_

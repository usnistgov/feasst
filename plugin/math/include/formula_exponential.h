
#ifndef FEASST_MATH_FORMULA_EXPONENTIAL_H_
#define FEASST_MATH_FORMULA_EXPONENTIAL_H_

#include <memory>
#include <string>
#include <sstream>
#include "utils/include/arguments.h"
#include "math/include/formula.h"

namespace feasst {

/**
  \f$ f(x) = f(x_0) \exp[A(x - x_0)^B] \f$
 */
class FormulaExponential : public Formula {
 public:
  /**
    args:
    - f0: function value at reference point (default: 0).
    - A: coefficient inside exponential (default: 1).
    - B: stretched exponential by power (default: 1).
   */
  FormulaExponential(const argtype& args = argtype());
  void set_f0(const double f0) { f0_ = f0; }
  void set_A(const double A) { A_ = A; }
  void set_B(const double B) { B_ = B; }
  double evaluate(const double x) const override;
  std::shared_ptr<Formula> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit FormulaExponential(std::istream& istr);
  virtual ~FormulaExponential() {}

 private:
  double f0_;
  double A_;
  double B_;
};

inline std::shared_ptr<FormulaExponential> MakeFormulaExponential(
    const argtype &args = argtype()) {
  return std::make_shared<FormulaExponential>(args);
}

}  // namespace feasst

#endif  // FEASST_MATH_FORMULA_EXPONENTIAL_H_


#ifndef FEASST_CORE_FORMULA_EXPONENTIAL_H_
#define FEASST_CORE_FORMULA_EXPONENTIAL_H_

#include <math.h>
#include "core/include/formula.h"

namespace feasst {

/**
  \f$ f(x) = f(x_0) \exp[A(x - x_0)^B] \f$
 */
class FormulaExponential : public Formula {
 public:
  FormulaExponential& set_f0(const double f0) { f0_ = f0; return *this; }
  FormulaExponential& set_x0(const double x0) { x0_ = x0; return *this; }
  FormulaExponential& set_A(const double A) { A_ = A; return *this; }
  FormulaExponential& set_B(const double B) { B_ = B; return *this; }

  double evaluate(const double x) {
    return f0_*exp(A_*pow(x - x0_, B_));
  }

  virtual ~FormulaExponential() {}

 private:
  double f0_;
  double x0_;
  double A_;
  double B_;
};

}  // namespace feasst

#endif  // FEASST_CORE_FORMULA_EXPONENTIAL_H_

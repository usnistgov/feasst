
#ifndef FEASST_CORE_FORMULA_POLYNOMIAL_H_
#define FEASST_CORE_FORMULA_POLYNOMIAL_H_

#include <vector>
#include <math.h>
#include "core/include/formula.h"

namespace feasst {

/**
  \f$ f(x) = \sum_{i=0}^n A_n(x - x_0)^i \f$
  where A_0 = f(x_0).
 */
class FormulaPolynomial : public Formula {
 public:
  FormulaPolynomial& set_A(const int index, const double A) {
    // resize A coefficients if necessary
    if (static_cast<int>(A_.size()) <= index) {
      A_.resize(index + 1);
    }
    A_[index] = A;
    return *this;
  }

  double evaluate(const double x) {
    double result = 0.;
    for (int i = 0; i < static_cast<int>(A_.size()); ++i) {
      result += A_[i]*pow(x - x0_, i);
    }
    return result;
  }

  virtual ~FormulaPolynomial() {}

 private:
  std::vector<double> A_;
};

}  // namespace feasst

#endif  // FEASST_CORE_FORMULA_POLYNOMIAL_H_

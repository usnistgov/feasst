
#ifndef FEASST_MATH_FORMULA_POLYNOMIAL_H_
#define FEASST_MATH_FORMULA_POLYNOMIAL_H_

#include <vector>
#include <math.h>
#include "math/include/formula.h"

namespace feasst {

/**
  \f$ f(x) = \sum_{i=0}^n A_n(x - x_0)^i \f$
  where A_0 = f(x_0).
 */
class FormulaPolynomial : public Formula {
 public:
  FormulaPolynomial(
    /**
     */
    const argtype& args = argtype()) : Formula(args) {}

  FormulaPolynomial& set_A(const int index, const double A) {
    // resize A coefficients if necessary
    if (static_cast<int>(A_.size()) <= index) {
      A_.resize(index + 1);
    }
    A_[index] = A;
    return *this;
  }

  double evaluate(const double x) override {
    double result = 0.;
    for (int i = 0; i < static_cast<int>(A_.size()); ++i) {
      result += A_[i]*pow(x - x0(), i);
    }
    return result;
  }

  std::shared_ptr<Formula> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  FormulaPolynomial(std::istream& istr);
  virtual ~FormulaPolynomial() {}

 private:
  const std::string class_name_ = "FormulaPolynomial";
  std::vector<double> A_;
};

inline std::shared_ptr<FormulaPolynomial> MakeFormulaPolynomial(
    const argtype &args = argtype()) {
  return std::make_shared<FormulaPolynomial>(args);
}

}  // namespace feasst

#endif  // FEASST_MATH_FORMULA_POLYNOMIAL_H_

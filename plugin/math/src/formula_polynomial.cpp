#include <cmath>
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/formula_polynomial.h"

namespace feasst {

FEASST_MAPPER(FormulaPolynomial,);

FormulaPolynomial::FormulaPolynomial(argtype args) : FormulaPolynomial(&args) {
  feasst_check_all_used(args); }
FormulaPolynomial::FormulaPolynomial(argtype * args) : Formula(args) {
  class_name_ = "FormulaPolynomial";
}

std::shared_ptr<Formula> FormulaPolynomial::create(std::istream& istr) const {
  return std::make_shared<FormulaPolynomial>(istr);
}

FormulaPolynomial::FormulaPolynomial(std::istream& istr)
  : Formula(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6937, "version mismatch: " << version);
  feasst_deserialize(&A_, istr);
}

void FormulaPolynomial::serialize(std::ostream& ostr) const {
  serialize_formula_(ostr);
  feasst_serialize_version(6937, ostr);
  feasst_serialize(A_, ostr);
}

double FormulaPolynomial::evaluate(const double x) const {
  double result = 0.;
  for (int i = 0; i < static_cast<int>(A_.size()); ++i) {
    result += A_[i]*std::pow(x - x0(), i);
  }
  return result;
}

double FormulaPolynomial::derivative(const double x) const {
  double result = 0.;
  for (int i = 1; i < static_cast<int>(A_.size()); ++i) {
    result += A_[i]*i*std::pow(x - x0(), i - 1);
  }
  return result;
}

}  // namespace feasst

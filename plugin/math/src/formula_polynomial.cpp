#include <cmath>
#include "math/include/formula_polynomial.h"
#include "utils/include/utils_io.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapFormulaPolynomial {
 public:
  MapFormulaPolynomial() {
    FormulaPolynomial().deserialize_map()["FormulaPolynomial"] =
      std::make_shared<FormulaPolynomial>();
  }
};

static MapFormulaPolynomial mapper_ = MapFormulaPolynomial();

std::shared_ptr<Formula> FormulaPolynomial::create(std::istream& istr) const {
  return std::make_shared<FormulaPolynomial>(istr);
}

FormulaPolynomial::FormulaPolynomial(std::istream& istr)
  : Formula(istr) {
  feasst_deserialize_version(istr);
  feasst_deserialize(&A_, istr);
}

void FormulaPolynomial::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_formula_(ostr);
  feasst_serialize_version(1, ostr);
  feasst_serialize(A_, ostr);
}

double FormulaPolynomial::evaluate(const double x) {
  double result = 0.;
  for (int i = 0; i < static_cast<int>(A_.size()); ++i) {
    result += A_[i]*std::pow(x - x0(), i);
  }
  return result;
}

}  // namespace feasst

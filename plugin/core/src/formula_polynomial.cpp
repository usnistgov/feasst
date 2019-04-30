
#include "core/include/formula_polynomial.h"
#include "core/include/utils_io.h"
#include "core/include/debug.h"

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

}  // namespace feasst

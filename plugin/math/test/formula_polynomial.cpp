#include <sstream>
#include "utils/test/utils.h"
#include "math/include/formula_polynomial.h"

namespace feasst {

TEST(FormulaPolynomial, serialize) {
  FormulaPolynomial formula;
  test_serialize<FormulaPolynomial, Formula>(formula, "FormulaPolynomial 5694 0 6937 0 ");
}

}  // namespace feasst

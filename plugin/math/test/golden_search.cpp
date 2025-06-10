#include "utils/test/utils.h"
#include "math/include/golden_search.h"
#include "math/include/formula_polynomial.h"

namespace feasst {

TEST(GoldenSearch, minimum) {
  // f(x) = (x - 2.0001)^2 + 1.2 = x^2 -4.0002x + 5.20040001
  auto poly = MakeFormulaPolynomial();
  poly->set_A(0, 5.20040001);
  poly->set_A(1, -4.0002);
  poly->set_A(2, 1.);
  const double res = GoldenSearch({
    {"lower", "-100"},
    {"upper", "100"},
    {"tolerance", str(1e-8)}}).minimum(poly.get());
  EXPECT_NEAR(res, 2.0001, 1e-6);
}

}  // namespace feasst

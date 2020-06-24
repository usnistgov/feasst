#include "utils/test/utils.h"
#include "math/include/solver_brent_dekker.h"
#include "math/include/formula_polynomial.h"

namespace feasst {

TEST(SolverBrentDekker, root) {
  auto poly = MakeFormulaPolynomial();
  poly->set_A(0, 2);
  poly->set_A(2, -1.);
  poly->set_A(3, 1.);
  const double res = SolverBrentDekker({
    {"lower", "-200"},
    {"upper", "300"},
    {"tolerance", "0.0025"}}).root(poly.get());
  EXPECT_NEAR(-1, res, 0.0025);
}

}  // namespace feasst

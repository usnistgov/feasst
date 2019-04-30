#include <sstream>
#include <gtest/gtest.h>
#include "core/include/formula_polynomial.h"

namespace feasst {

TEST(FormulaPolynomial, serialize) {
  FormulaPolynomial formula;
  std::stringstream ss, ss2;
  formula.serialize(ss);
  EXPECT_EQ("FormulaPolynomial 1 0 1 0 ", ss.str());
  std::shared_ptr<Formula> formula2 = FormulaPolynomial().deserialize(ss);
  formula2->serialize(ss2);
  EXPECT_EQ(ss2.str(), ss.str());
}

}  // namespace feasst

#include <cmath>
#include "utils/test/utils.h"
#include "flat_histogram/include/ln_probability.h"

namespace feasst {

TEST(LnProbability, serialize) {
  LnProbability ln_prob;
  ln_prob.resize(5);
  ln_prob.add(0, 1.);
  ln_prob.add(1, 0.5);
  EXPECT_DOUBLE_EQ(ln_prob.value(0), 1);
  EXPECT_DOUBLE_EQ(ln_prob.value(1), 0.5);
  EXPECT_DOUBLE_EQ(ln_prob.delta(1), -0.5);
  EXPECT_DOUBLE_EQ(ln_prob.delta_values()[1], -0.5);
  EXPECT_TRUE(std::isnan(ln_prob.delta(0)));
  LnProbability ln_prob2 = test_serialize(ln_prob, "885 5 1 0.5 0 0 0 ");
}

}  // namespace feasst

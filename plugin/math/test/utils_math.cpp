#include <vector>
#include <gtest/gtest.h>
#include "math/include/utils_math.h"

namespace feasst {

TEST(UtilsMath, sgn) {
  EXPECT_EQ(sgn(-12.5), -1);
  EXPECT_EQ(sgn(5./2.), 1);
}

TEST(UtilsMath, cumulative_probability) {
  std::vector<double> weights = {1, 2, 7};
  std::vector<double> cpdf = cumulative_probability(weights);
  EXPECT_NEAR(0.1, cpdf[0], NEAR_ZERO);
  EXPECT_NEAR(0.1 + 0.2, cpdf[1], NEAR_ZERO);
  EXPECT_NEAR(0.1 + 0.2 + 0.7, cpdf[2], NEAR_ZERO);
}

}  // namespace feasst

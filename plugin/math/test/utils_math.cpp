#include <vector>
#include <gtest/gtest.h>
#include "math/include/utils_math.h"
#include "utils/include/utils_io.h"

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

TEST(UtilsMath, minimum) {
  std::vector<double> data = {10, 9, 8, 8.5, 7, 0.122, 6, 8, 10};
  EXPECT_EQ(0.122, minimum(data));
  std::vector<int> mins = local_minimum_indices(data, 1);
  EXPECT_EQ(2, mins.size());
  EXPECT_EQ(2, mins[0]);
  EXPECT_EQ(5, mins[1]);
  for (int smooth : {3, 11}) {
    mins = local_minimum_indices(data, smooth);
    EXPECT_EQ(1, mins.size());
    EXPECT_EQ(5, mins[0]);
  }
}

}  // namespace feasst

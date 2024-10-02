#include <vector>
#include "utils/test/utils.h"
#include "utils/include/io.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"

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

TEST(UtilsMath, spherical_shell_volume) {
  EXPECT_NEAR(spherical_shell_volume(1.7, 3.2, 2), 23.090706003884986, NEAR_ZERO);
  EXPECT_NEAR(spherical_shell_volume(1.7, 3.2, 3), 116.67875115432494, NEAR_ZERO);
}

TEST(UtilsMath, is_in_interval) {
  EXPECT_TRUE(is_in_interval(5, 6, 4));
  EXPECT_TRUE(is_in_interval(5, 4.1, 6.1));
  EXPECT_FALSE(is_in_interval(5, 7.1, 6.1));
  EXPECT_FALSE(is_in_interval(5, 6.1, 7.1));
}

TEST(UtilsMath, factorial) {
  EXPECT_EQ(1, factorial(0));
  EXPECT_EQ(1, factorial(1));
  EXPECT_EQ(2, factorial(2));
  EXPECT_EQ(6, factorial(3));
  EXPECT_EQ(24, factorial(4));
  EXPECT_NEAR(3628800, factorial(10), 1e-9);
  TRY(
    factorial(-1);
    CATCH_PHRASE("not implemented");
  );
}

}  // namespace feasst

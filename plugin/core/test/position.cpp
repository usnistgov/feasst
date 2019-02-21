#include <gtest/gtest.h>
#include <math.h>
#include "core/test/position_test.h"
#include "core/include/constants.h"
#include "core/include/debug.h"

namespace feasst {

TEST(Position, getset) {
  Position pos = default_position();
  EXPECT_EQ(3, pos.size());
  std::vector<double> x = {3.5, 796.4, -45.4};
  pos.set_vector(x);
  std::vector<double> x2 = pos.coord();
  EXPECT_EQ(x, x2);
  EXPECT_EQ(3, x.size());
  EXPECT_EQ(pos.coord(1), 796.4);
  EXPECT_NEAR(pos.dot_product(pos), 636326.37, NEAR_ZERO);
  EXPECT_NEAR(pos.dot_product(pos), pos.squared_distance(), NEAR_ZERO);
}

}  // namespace feasst

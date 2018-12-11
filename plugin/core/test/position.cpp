#include <gtest/gtest.h>
#include "core/include/position.h"

TEST(Position, getset) {
  feasst::Position pos;
  pos.set_to_origin_3D();
  EXPECT_EQ(3, pos.size());
  std::vector<double> x = {3.5, 796.4, -45.4};
  pos.set_vector(x);
  std::vector<double> x2 = pos.coord();
  EXPECT_EQ(x, x2);
  EXPECT_EQ(3, x.size());
  EXPECT_EQ(pos.coord(1), 796.4);
  EXPECT_NEAR(pos.dot_product(pos), 636326.37, 1e-15);
  EXPECT_NEAR(pos.dot_product(pos), pos.squared_distance(), 1e-15);
}

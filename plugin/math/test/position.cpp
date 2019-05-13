#include "utils/test/utils.h"
#include <math.h>
#include "system/test/position_test.h"
#include "math/include/constants.h"
#include "utils/include/debug.h"

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

  // serialize
  Position pos2 = test_serialize(pos);
  EXPECT_EQ(pos.coord(), pos2.coord());
}

TEST(Position, cross_product) {
  auto a = Position().set_vector({1, 2, 3});
  auto b = Position().set_vector({3.5, 3.5, 3.5});
  Position c = a.cross_product(b);
  EXPECT_NEAR(-3.5, c.coord(0), 1e-10);
  EXPECT_NEAR(7., c.coord(1), 1e-10);
  EXPECT_NEAR(-3.5, c.coord(2), 1e-10);
}

TEST(Position, nearest_distance_to_axis) {
  auto x0 = Position().set_vector({1., 2., 3.});
  auto x1 = Position().set_vector({2., 3., 4.});
  auto x2 = Position().set_vector({3., 4., 5.});
  EXPECT_NEAR(x0.nearest_distance_to_axis(x1, x2), 0, NEAR_ZERO);
}

TEST(Position, nearest_distance_to_axis2) {
  auto x0 = Position().set_vector({2.354, -1, 1.});
  auto x1 = Position().set_vector({0, -1., 0.});
  auto x2 = Position().set_vector({0., -1., 1});
  EXPECT_NEAR(x0.nearest_distance_to_axis(x1, x2), 2.354, NEAR_ZERO);
}

}  // namespace feasst

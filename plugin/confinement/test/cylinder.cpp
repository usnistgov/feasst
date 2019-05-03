#include <gtest/gtest.h>
#include "confinement/include/cylinder.h"
#include "core/include/debug.h"

namespace feasst {

TEST(Shape, Cylinder) {
  Cylinder cylinder(
    {{"radius", "2"}},
    Position().set_vector({0, 0, 0}),
    Position().set_vector({0, 0, 1})
  );
  Position point;
  point.set_vector({1.5, 0, 13.535});
  EXPECT_NEAR(-0.5, cylinder.nearest_distance(point), 1e-15);
  EXPECT_TRUE(cylinder.is_inside(point));
  EXPECT_TRUE(cylinder.is_inside(point, 0.9999));
  EXPECT_FALSE(cylinder.is_inside(point, 1.00001));

  std::stringstream ss, ss2;
  cylinder.serialize(ss);
  std::shared_ptr<Shape> cylinder2 = cylinder.deserialize(ss);
  EXPECT_TRUE(cylinder2->is_inside(point, 0.9999));
  EXPECT_FALSE(cylinder2->is_inside(point, 1.00001));
  cylinder2->serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst

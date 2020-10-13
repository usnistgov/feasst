#include "utils/test/utils.h"
#include "utils/include/debug.h"
#include "shape/include/half_space_tilted.h"

namespace feasst {

TEST(HalfSpaceTilted, serialize) {
  HalfSpaceTilted half_space(Position({0, 0, 1}), Position({0, 0, 12}));
  EXPECT_NEAR(half_space.unit_normal().coord(0), 0, NEAR_ZERO);
  EXPECT_NEAR(half_space.unit_normal().coord(1), 0, NEAR_ZERO);
  EXPECT_NEAR(half_space.unit_normal().coord(2), 1, NEAR_ZERO);
  EXPECT_NEAR(half_space.distance_from_origin(), 1, NEAR_ZERO);

  Position point;
  point.set_vector({15, -56.54, 2.});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({15, -56.54, 1.0000001});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({15, -56.54, 0.9999999});
  EXPECT_FALSE(half_space.is_inside(point));
  point.set_vector({15, -56.54, 1.5000001});
  EXPECT_TRUE(half_space.is_inside(point, 1.));
  point.set_vector({15, -56.54, 1.4999999});
  EXPECT_FALSE(half_space.is_inside(point, 1.));

  // serialize
  auto half_space3 = test_serialize<HalfSpaceTilted, Shape>(half_space);
  point.set_vector({15, -56.54, 1.4999999999});
  EXPECT_FALSE(half_space3->is_inside(point, 1.));

  HalfSpaceTilted half_space2(Position({0, 0, 3}), Position({0, 0, -12}));
  EXPECT_NEAR(half_space2.unit_normal().coord(0), 0, NEAR_ZERO);
  EXPECT_NEAR(half_space2.unit_normal().coord(1), 0, NEAR_ZERO);
  EXPECT_NEAR(half_space2.unit_normal().coord(2), -1, NEAR_ZERO);
  EXPECT_NEAR(half_space2.distance_from_origin(), -3, NEAR_ZERO);
  EXPECT_TRUE(half_space2.is_inside(point));
  point.set_vector({15, -56.54, 3.000000000000001});
  EXPECT_FALSE(half_space2.is_inside(point));
  point.set_vector({15, -56.54, 2.9999999999});
  EXPECT_TRUE(half_space2.is_inside(point));
  point.set_vector({15, -56.54, -50});
  EXPECT_TRUE(half_space2.is_inside(point));

  HalfSpaceTilted half_space4(Position({0, 0, 3}), Position({0, 0, 0}));
  EXPECT_EQ(half_space4.unit_normal().coord(2),
            half_space2.unit_normal().coord(2));
  EXPECT_EQ(half_space4.distance_from_origin(),
            half_space2.distance_from_origin());
}

}  // namespace feasst

#include "utils/test/utils.h"
#include "utils/include/debug.h"
#include "shape/include/cuboid.h"

namespace feasst {

TEST(Cuboid, volume) {
  auto cuboid = MakeCuboid({{"side_lengths", "1,2,3"}, {"center", "0,0,6"}});
  Cuboid cuboid2 = test_serialize(*cuboid);
  EXPECT_NEAR(cuboid2.volume(), 6, NEAR_ZERO);
  EXPECT_NEAR(cuboid2.surface_area(), 22, NEAR_ZERO);
  EXPECT_EQ(cuboid2.center().coord(0), 0);
  EXPECT_EQ(cuboid2.center().coord(1), 0);
  EXPECT_EQ(cuboid2.center().coord(2), 6);
  auto cube = MakeCuboid({{"cubic_side_length", "3"}});
  EXPECT_NEAR(cube->volume(), 27, NEAR_ZERO);
  EXPECT_NEAR(cube->surface_area(), 54, NEAR_ZERO);
}

}  // namespace feasst

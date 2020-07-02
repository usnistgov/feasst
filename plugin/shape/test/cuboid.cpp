#include "utils/test/utils.h"
#include "utils/include/debug.h"
#include "shape/include/cuboid.h"

namespace feasst {

TEST(Cuboid, volume) {
  Cuboid cuboid(Position({1, 2, 3}), Position({0, 0, 0}));
  Cuboid cuboid2 = test_serialize(cuboid);
  EXPECT_NEAR(cuboid2.volume(), 6, NEAR_ZERO);
  EXPECT_NEAR(cuboid2.surface_area(), 22, NEAR_ZERO);
  auto cube = MakeCube(3);
  EXPECT_NEAR(cube->volume(), 27, NEAR_ZERO);
  EXPECT_NEAR(cube->surface_area(), 54, NEAR_ZERO);
}

}  // namespace feasst

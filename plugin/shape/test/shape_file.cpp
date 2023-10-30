#include "utils/test/utils.h"
#include "utils/include/debug.h"
#include "shape/include/shape_file.h"

namespace feasst {

TEST(ShapeFile, shape) {
  auto shape = MakeShapeFile({{"shape_file", "../plugin/shape/test/data/shape.txt"}});
  ShapeFile shape2 = test_serialize(*shape);
//  EXPECT_NEAR(cuboid2.volume(), 6, NEAR_ZERO);
//  EXPECT_NEAR(cuboid2.surface_area(), 22, NEAR_ZERO);
//  EXPECT_EQ(cuboid2.center().coord(0), 0);
//  EXPECT_EQ(cuboid2.center().coord(1), 0);
//  EXPECT_EQ(cuboid2.center().coord(2), 6);
//  auto cube = MakeCuboid({{"cubic_side_length", "3"}});
//  EXPECT_NEAR(cube->volume(), 27, NEAR_ZERO);
//  EXPECT_NEAR(cube->surface_area(), 54, NEAR_ZERO);
}

}  // namespace feasst

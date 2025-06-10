#include "utils/test/utils.h"
#include "utils/include/debug.h"
#include "math/include/accumulator.h"
#include "shape/include/supertoroid.h"

namespace feasst {

TEST(Shape, Supertoroid) {
  auto shape = MakeSupertoroid({{"a1", "2"}, {"center", "0,0,0"}});
  std::shared_ptr<Shape> shape2 = test_serialize<Supertoroid, Shape>(*shape);
  Position point;
  point.set_vector({1.5, 0, 0});
  EXPECT_TRUE(shape2->is_inside(point, 0.9999));
  point.set_coord(0, 4);
  EXPECT_FALSE(shape2->is_inside(point, 0.9999));
}

}  // namespace feasst

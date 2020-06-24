#include "utils/test/utils.h"
#include "shape/include/slab.h"
#include "utils/include/debug.h"

namespace feasst {

TEST(Slab, serialize) {
  Slab slab({
    {"dimension", "2"},
    {"bound0", "3"},
    {"bound1", "1"},
  });
  Position point;
  point.set_vector({15, -56.54, 1.5});
  EXPECT_NEAR(-0.5, slab.nearest_distance(point), 1e-15);
  EXPECT_TRUE(slab.is_inside(point));
  EXPECT_TRUE(slab.is_inside(point, 0.9999));
  EXPECT_FALSE(slab.is_inside(point, 1.00001));

  std::shared_ptr<Shape> slab2 = test_serialize<Slab, Shape>(slab);
  EXPECT_TRUE(slab2->is_inside(point, 0.9999));
  EXPECT_FALSE(slab2->is_inside(point, 1.00001));
}

}  // namespace feasst

#include <gtest/gtest.h>
#include "confinement/include/slab.h"
#include "core/include/debug.h"

namespace feasst {

TEST(Shape, Slab) {
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

  std::stringstream ss, ss2;
  slab.serialize(ss);
  std::shared_ptr<Shape> slab2 = slab.deserialize(ss);
  EXPECT_TRUE(slab2->is_inside(point, 0.9999));
  EXPECT_FALSE(slab2->is_inside(point, 1.00001));
  slab2->serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst

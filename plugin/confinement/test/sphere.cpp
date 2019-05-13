#include "utils/test/utils.h"
#include "confinement/include/sphere.h"
#include "utils/include/debug.h"

namespace feasst {

TEST(Shape, Sphere) {
  Sphere sphere(
    {{"radius", "2"}},
    Position().set_vector({0, 0, 0})
  );
  Position point;
  point.set_vector({1.5, 0, 0});
  EXPECT_NEAR(-0.5, sphere.nearest_distance(point), 1e-15);
  EXPECT_TRUE(sphere.is_inside(point));
  EXPECT_TRUE(sphere.is_inside(point, 0.9999));
  EXPECT_FALSE(sphere.is_inside(point, 1.00001));

  std::shared_ptr<Shape> sphere2 = test_serialize<Sphere, Shape>(sphere);
  EXPECT_TRUE(sphere2->is_inside(point, 0.9999));
  EXPECT_FALSE(sphere2->is_inside(point, 1.00001));
}

}  // namespace feasst

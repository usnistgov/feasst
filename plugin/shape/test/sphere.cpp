#include "utils/test/utils.h"
#include "utils/include/debug.h"
//#include "math/include/histogram.h"
#include "math/include/accumulator.h"
#include "shape/include/sphere.h"

namespace feasst {

TEST(Shape, Sphere) {
  auto sphere = MakeSphere({{"radius", "2"}, {"center", "0,0,0"}});
  Position point;
  point.set_vector({1.5, 0, 0});
  EXPECT_NEAR(-0.5, sphere->nearest_distance(point), 1e-15);
  EXPECT_TRUE(sphere->is_inside(point));
  EXPECT_TRUE(sphere->is_inside(point, 0.9999));
  EXPECT_FALSE(sphere->is_inside(point, 1.00001));

  std::shared_ptr<Shape> sphere2 = test_serialize<Sphere, Shape>(*sphere);
  EXPECT_TRUE(sphere2->is_inside(point, 0.9999));
  EXPECT_FALSE(sphere2->is_inside(point, 1.00001));

  auto sphere3 = MakeSphere({{"radius", "2"}});
  EXPECT_EQ(sphere3->center().size(), 3);
  EXPECT_EQ(sphere3->center().coord(0), 0);
}

TEST(Sphere, surface_mesh) {
  std::vector<Position> mesh;
  MakeSphere({{"radius", "3"}})->surface_mesh(1000, &mesh);
  Accumulator av_x, av_y, av_z;
  for (const Position& pos1 : mesh) {
    EXPECT_NEAR(pos1.distance(), 3, NEAR_ZERO);
    av_x.accumulate(pos1.coord(0));
    av_y.accumulate(pos1.coord(1));
    av_z.accumulate(pos1.coord(2));
  }
  EXPECT_LT(av_x.average(), 0.0003);
  EXPECT_LT(av_y.average(), 0.0003);
  EXPECT_LT(av_z.average(), 0.0003);
}

}  // namespace feasst

#include <gtest/gtest.h>
#include "confinement/include/shape.h"
#include "core/include/debug.h"

TEST(Shape, HalfSpace) {
  auto half_space = feasst::HalfSpace()
    .set_dimension(2)
    .set_intersection(1)
    .set_direction(1);
  try {
    half_space.set_direction(0.);
    CATCH_PHRASE("direction cannot be infinitesimal");
  }
  feasst::Position point;
  point.set_vector({15, -56.54, 2.});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({15, -56.54, 1.000000000000001});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({15, -56.54, 0.9999999});
  EXPECT_FALSE(half_space.is_inside(point));
  point.set_vector({15, -56.54, 1.500000000001});
  EXPECT_TRUE(half_space.is_inside(point, 1.));
  point.set_vector({15, -56.54, 1.4999999999});
  EXPECT_FALSE(half_space.is_inside(point, 1.));

  auto half_space2 = half_space;
  half_space2.set_intersection(3).set_direction(-1);
  EXPECT_TRUE(half_space2.is_inside(point));
  point.set_vector({15, -56.54, 3.000000000000001});
  EXPECT_FALSE(half_space2.is_inside(point));
  point.set_vector({15, -56.54, 2.9999999999});
  EXPECT_TRUE(half_space2.is_inside(point));

  feasst::ShapeIntersect slit(std::make_shared<feasst::HalfSpace>(half_space),
                              std::make_shared<feasst::HalfSpace>(half_space2));
  point.set_vector({15, -56.54, 1.5});
  EXPECT_NEAR(-0.5, slit.nearest_distance(point), 1e-15);
  EXPECT_TRUE(slit.is_inside(point));
  EXPECT_TRUE(slit.is_inside(point, 0.9999));
  EXPECT_FALSE(slit.is_inside(point, 1.00001));
}

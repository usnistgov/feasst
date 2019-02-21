#include <gtest/gtest.h>
#include <math.h>
#include "core/include/matrix.h"
#include "core/include/constants.h"
#include "core/include/debug.h"

namespace feasst {

TEST(Matrix, axis_angle) {
  RotationMatrix mat;
  Position axis;
  axis.set_vector({2., 0., 0.});
  mat.axis_angle(axis, 90);
  EXPECT_NEAR(mat.determinant(), 1., NEAR_ZERO);
  EXPECT_EQ(mat.multiply(axis).coord(), axis.coord());
  Position point1;
  point1.set_vector({0., 2., 0.});
  Position point2;
  point2.set_vector({0., 0., 2.});
  EXPECT_TRUE(point2.is_equal(mat.multiply(point1)));
  mat.axis_angle(axis, -135);
  point2.set_vector({0., -sqrt(2.), -sqrt(2.)});
  EXPECT_TRUE(point2.is_equal(mat.multiply(point1)));
}

}  // namespace feasst

#include <cmath>
#include "utils/test/utils.h"
#include "math/include/matrix.h"
#include "math/include/constants.h"
#include "utils/include/debug.h"

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
  point2.set_vector({0., -std::sqrt(2.), -std::sqrt(2.)});
  EXPECT_TRUE(point2.is_equal(mat.multiply(point1)));

  // For a rotation matrix, the inverse should be equal to the transpose
  RotationMatrix mat_t = mat, mat_inv = mat;
  mat_t.transpose();
  mat_inv.invert();
  EXPECT_TRUE(mat_t.is_equal(mat_inv));
  EXPECT_FALSE(mat.is_equal(mat_inv));
}

}  // namespace feasst

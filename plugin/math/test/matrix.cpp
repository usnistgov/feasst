#include <cmath>
#include "utils/test/utils.h"
#include "utils/include/debug.h"
#include "math/include/matrix.h"
#include "math/include/constants.h"

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

TEST(Matrix, 2d) {
  RotationMatrix mat;
  Position axis;
  axis.set_vector({0., 0.});
  mat.axis_angle(axis, 90);
  EXPECT_NEAR(mat.determinant(), 1., NEAR_ZERO);
  Position point1;
  point1.set_vector({2., 0.});
  Position point2;
  point2.set_vector({0., 2.});
  EXPECT_TRUE(point2.is_equal(mat.multiply(point1)));
  mat.axis_angle(axis, -135);
  point2.set_vector({-std::sqrt(2.), -std::sqrt(2.)});
  EXPECT_TRUE(point2.is_equal(mat.multiply(point1)));

  // For a rotation matrix, the inverse should be equal to the transpose
  RotationMatrix mat_t = mat, mat_inv = mat;
  mat_t.transpose();
  mat_inv.invert();
  EXPECT_TRUE(mat_t.is_equal(mat_inv));
  EXPECT_FALSE(mat.is_equal(mat_inv));
}

TEST(Matrix, multiply) {
  Matrix A({{1, 2, 3}, {4, 5, 6}});
  Matrix B({{7, 8}, {9, 10}, {11, 12}});
  Matrix C = A.multiply(B);
  EXPECT_EQ(C.value(0, 0), 58);
  EXPECT_EQ(C.value(0, 1), 64);
  EXPECT_EQ(C.value(1, 0), 139);
  EXPECT_EQ(C.value(1, 1), 154);
}

TEST(Matrix, is_identity) {
  Matrix mat({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
  EXPECT_TRUE(mat.is_identity());
}

TEST(Matrix, skew_symmetric_cross_product) {
  Position a({{"x", "1.6"}, {"y", "2.5"}, {"z", "3.4"}});
  Position b({{"x", "4.3"}, {"y", "5.2"}, {"z", "6.1"}});
  Position c = a.cross_product(b);
  Matrix ssc;
  ssc.skew_symmetric_cross_product(a);
  Position c2 = ssc.multiply(b);
  EXPECT_TRUE(c.is_equal(c2, NEAR_ZERO));
}

TEST(Matrix, identiy) {
  Matrix mat;
  mat.identity(3);
  EXPECT_NEAR(mat.value(0, 0), 1., NEAR_ZERO);
  EXPECT_NEAR(mat.value(0, 1), 0., NEAR_ZERO);
  EXPECT_NEAR(mat.value(0, 2), 0., NEAR_ZERO);
  EXPECT_NEAR(mat.value(1, 0), 0., NEAR_ZERO);
  EXPECT_NEAR(mat.value(1, 1), 1., NEAR_ZERO);
  EXPECT_NEAR(mat.value(1, 2), 0., NEAR_ZERO);
  EXPECT_NEAR(mat.value(2, 0), 0., NEAR_ZERO);
  EXPECT_NEAR(mat.value(2, 1), 0., NEAR_ZERO);
  EXPECT_NEAR(mat.value(2, 2), 1., NEAR_ZERO);
}

TEST(Matrix, add) {
  Matrix a;
  a.identity(2);
  a.add(a);
  EXPECT_NEAR(a.value(0, 0), 2., NEAR_ZERO);
  EXPECT_NEAR(a.value(0, 1), 0., NEAR_ZERO);
  EXPECT_NEAR(a.value(1, 0), 0., NEAR_ZERO);
  EXPECT_NEAR(a.value(1, 1), 2., NEAR_ZERO);
}

TEST(Matrix, vector_onto_vector) {
  Position a({{"x", "1.6"}, {"y", "2.5"}, {"z", "3.4"}});
  Position b({{"x", "4.3"}, {"y", "5.2"}, {"z", "6.1"}});
  RotationMatrix rot;
  rot.vector_onto_vector(a, b);
  Position c = rot.multiply(a);
  c.multiply(b.distance()/a.distance());
  EXPECT_TRUE(b.is_equal(c, NEAR_ZERO));
}

}  // namespace feasst

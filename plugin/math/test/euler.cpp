#include "utils/test/utils.h"
#include "math/include/euler.h"
#include "math/include/random_mt19937.h"

namespace feasst {

TEST(Euler, compute_rotation_matrix) {
  RotationMatrix mat;
  //mat.set_size(3, 3);
  Position x({0.6, 0.6, 0.6});
  Euler euler(PI/2., 0., 0.);
  euler.compute_rotation_matrix(&mat);
  DEBUG(mat.str());
  Position x_n = mat.multiply(x);
  EXPECT_NEAR(x_n.coord(0), -0.6, NEAR_ZERO);
  EXPECT_NEAR(x_n.coord(1), 0.6, NEAR_ZERO);
  EXPECT_NEAR(x_n.coord(2), 0.6, NEAR_ZERO);
  euler.set(0., PI/2., 0.);
  euler.compute_rotation_matrix(&mat);
  DEBUG(mat.str());
  x_n = mat.multiply(x);
  EXPECT_NEAR(x_n.coord(0), 0.6, NEAR_ZERO);
  EXPECT_NEAR(x_n.coord(1), -0.6, NEAR_ZERO);
  EXPECT_NEAR(x_n.coord(2), 0.6, NEAR_ZERO);
  euler.set(0., 0., PI/2.);
  euler.compute_rotation_matrix(&mat);
  DEBUG(mat.str());
  x_n = mat.multiply(x);
  EXPECT_NEAR(x_n.coord(0), -0.6, NEAR_ZERO);
  EXPECT_NEAR(x_n.coord(1), 0.6, NEAR_ZERO);
  EXPECT_NEAR(x_n.coord(2), 0.6, NEAR_ZERO);
  euler.set(PI, 0., 0);
  euler.compute_rotation_matrix(&mat);
  DEBUG(mat.str());
  x_n = mat.multiply(x);
  EXPECT_NEAR(x_n.coord(0), -0.6, NEAR_ZERO);
  EXPECT_NEAR(x_n.coord(1), -0.6, NEAR_ZERO);
  EXPECT_NEAR(x_n.coord(2), 0.6, NEAR_ZERO);
  euler.set(30./180.*PI, 30./180.*PI, 30./180.*PI);
  euler.compute_rotation_matrix(&mat);
  DEBUG(mat.str());
  x_n = mat.multiply(x);
  EXPECT_NEAR(x_n.coord(0), -0.01471143, 1e-8);
  EXPECT_NEAR(x_n.coord(1), 0.46471143, 1e-8);
  EXPECT_NEAR(x_n.coord(2), 0.92942286, 1e-8);

  // check inverse rotation
  Matrix transpose = mat;
  transpose.transpose();
  Matrix identity = mat.multiply(transpose);
  EXPECT_TRUE(identity.is_identity());

  // obtain Eulers from rotation matrix and vice versa
  RandomMT19937 random;
  Euler euler1, euler2;
  RotationMatrix rot;
  Position pos({0, 0, 1});
  for (int i = 0; i < 10; ++i) {
  //for (int i = 0; i < 1000; ++i) {
    euler1.set(PI*(2*random.uniform() - 1),
               PI*(random.uniform()),
               PI*(2*random.uniform() - 1));
    euler1.compute_rotation_matrix(&rot);
    euler2.set(rot);
    EXPECT_TRUE(euler1.is_equal(euler2));
    if (!euler1.is_equal(euler2)) {
      INFO(euler1.str());
      INFO(euler2.str());
    }
    Position npos = rot.multiply(pos);
    //std::cout << npos.coord(0) << " " << npos.coord(1) << " " << npos.coord(2) << std::endl;
  }
}

}  // namespace feasst

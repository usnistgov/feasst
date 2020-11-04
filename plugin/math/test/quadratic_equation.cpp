#include "utils/test/utils.h"
#include "math/include/quadratic_equation.h"

namespace feasst {

TEST(Math, quadratic_equation) {
  double discriminant = -1, root1 = -1, root2 = -1;

  quadratic_equation(0., 2., 1., &discriminant, &root1, &root2);
  EXPECT_NEAR(4, discriminant, NEAR_ZERO);
  EXPECT_NEAR(-0.5, root1, NEAR_ZERO);
  EXPECT_NEAR(-1, root2, NEAR_ZERO);

  quadratic_equation(1., 0., -4., &discriminant, &root1, &root2);
  EXPECT_NEAR(16., discriminant, NEAR_ZERO);
  EXPECT_NEAR(2, root1, NEAR_ZERO);
  EXPECT_NEAR(-2, root2, NEAR_ZERO);

  quadratic_equation(1., 5., 6., &discriminant, &root1, &root2);
  EXPECT_NEAR(1., discriminant, NEAR_ZERO);
  EXPECT_NEAR(-2, root1, NEAR_ZERO);
  EXPECT_NEAR(-3, root2, NEAR_ZERO);

  quadratic_equation(1., 3., 3., &discriminant, &root1, &root2);
  EXPECT_NEAR(-3, discriminant, NEAR_ZERO);
  EXPECT_NEAR(-2, root1, NEAR_ZERO);
  EXPECT_NEAR(-3, root2, NEAR_ZERO);
}

}  // namespace feasst

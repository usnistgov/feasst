#include "utils/test/utils.h"
#include "shape/include/formula_sine_wave.h"

namespace feasst {

TEST(FormulaSineWave, nearest) {
  FormulaSineWave sine;
  EXPECT_NEAR(sine.nearest_minimum(0), -PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_maximum(0), PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_minimum(PI/2 - 0.01), -PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_minimum(PI/2 + 0.01), 3*PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_maximum(3*PI/2 - 0.01), PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_maximum(3*PI/2 + 0.01), 5*PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_minimum(5*PI/2 - 0.01), 3*PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_minimum(5*PI/2 + 0.01), 7*PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_maximum(7*PI/2 - 0.01), 5*PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_maximum(7*PI/2 + 0.01), 9*PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_minimum(-3*PI/2 + 0.01), -PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_minimum(-3*PI/2 - 0.01), -5*PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_maximum(-PI/2 + 0.01), PI/2, NEAR_ZERO);
  EXPECT_NEAR(sine.nearest_maximum(-PI/2 - 0.01), -3*PI/2, NEAR_ZERO);

  FormulaSineWave sine2({{"width", "1"}, {"phase", "0.1"}});
  EXPECT_NEAR(sine2.nearest_minimum(0), -1./4 + 0.1, NEAR_ZERO);
  EXPECT_NEAR(sine2.nearest_maximum(0), 1./4 + 0.1, NEAR_ZERO);
  EXPECT_NEAR(sine2.nearest_minimum(1./4 + 0.1 - 0.01), -1./4 + 0.1, NEAR_ZERO);
  EXPECT_NEAR(sine2.nearest_minimum(1./4 + 0.1 + 0.01), 3./4 + 0.1, NEAR_ZERO);
  EXPECT_NEAR(sine2.nearest_maximum(3./4 + 0.1 - 0.01), 1./4 + 0.1, NEAR_ZERO);
  EXPECT_NEAR(sine2.nearest_maximum(3./4 + 0.1 + 0.01), 5./4 + 0.1, NEAR_ZERO);
  EXPECT_NEAR(sine2.nearest_minimum(5./4 + 0.1 - 0.01), 3./4 + 0.1, NEAR_ZERO);
  EXPECT_NEAR(sine2.nearest_minimum(5./4 + 0.1 + 0.01), 7./4 + 0.1, NEAR_ZERO);
  EXPECT_NEAR(sine2.nearest_minimum(-3./4 + 0.1 + 0.01), -1./4 + 0.1, NEAR_ZERO);
  EXPECT_NEAR(sine2.nearest_minimum(-3./4 + 0.1 - 0.01), -5./4 + 0.1, NEAR_ZERO);

  FormulaSineWave sine3({{"width", "8"}, {"phase", "4"}});
  EXPECT_NEAR(sine3.nearest_minimum(0), 2, NEAR_ZERO);
  EXPECT_NEAR(sine3.nearest_maximum(0), -2, NEAR_ZERO);
  EXPECT_NEAR(sine3.nearest_minimum(5.999), 2, NEAR_ZERO);
  EXPECT_NEAR(sine3.nearest_minimum(6.001), 10, NEAR_ZERO);
  EXPECT_NEAR(sine3.nearest_maximum(1.999), -2, NEAR_ZERO);
  EXPECT_NEAR(sine3.nearest_maximum(2.001), 6, NEAR_ZERO);
}

}  // namespace feasst

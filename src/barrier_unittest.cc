#include <gtest/gtest.h>
#include "barrier.h"

// create barrier walls for x >= 5 and x <= -5 
TEST(Barrier, planar) {
  feasst::Barrier barrier;
  barrier.addOrthogonalPlanar(5, 1, 0);
  barrier.addOrthogonalPlanar(-5, -1, 0);
  
  vector<double> x(3, 0.);
  EXPECT_NEAR(0., barrier.potential(x, 1.), feasst::DTOL);
  
  x[0] = 4.499;
  EXPECT_NEAR(0., barrier.potential(x, 1.), feasst::DTOL);
  
  x[0] = 4.501;
  EXPECT_LT(1e100, barrier.potential(x, 1.));

  x[0] = -4.499;
  EXPECT_NEAR(0., barrier.potential(x, 1.), feasst::DTOL);
  
  x[0] = -4.501;
  EXPECT_LT(1e100, barrier.potential(x, 1.));

}


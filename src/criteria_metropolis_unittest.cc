#include <gtest/gtest.h>
#include <limits.h>
#include "criteria.h"
#include "criteria_metropolis.h"

TEST(Criteria, Metropolis) {
  const double beta = 1e-15, activ = 1;
  feasst::CriteriaMetropolis c(beta, activ);
  //feasst::CriteriaMetropolis c(beta, activ);
  EXPECT_EQ(1, c.accept(0., 0., "move", 0));
}
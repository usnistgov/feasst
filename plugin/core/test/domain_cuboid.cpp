#include <gtest/gtest.h>
#include "core/include/domain_cuboid.h"
#include "core/include/constants.h"

TEST(DomainCuboid, min_side_length) {
  feasst::DomainCuboid domain;
  domain.set_cubic(10);
  EXPECT_NEAR(10, domain.min_side_length(), 1e-15);
  domain.set_side_length(1, 5.234);
  EXPECT_NEAR(5.234, domain.min_side_length(), 1e-15);
}

// compute the volume of a sphere
TEST(DomainCuboid, random_position) {
  feasst::DomainCuboid domain;
  feasst::Random random;
  domain.set_cubic(1);
  int inside = 0;
  const int trials = 1e3;
  for (int trial = 0; trial < trials; ++trial) {
    if (domain.random_position(&random).squared_distance() <= 1) {
      ++inside;
    }
  }
  EXPECT_NEAR(static_cast<double>(inside)/trials, feasst::PI/6, 0.06);
}

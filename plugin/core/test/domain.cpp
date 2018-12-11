#include <gtest/gtest.h>
#include "core/include/domain_cuboid.h"

TEST(DomainCuboid, Volume) {
  feasst::DomainCuboid domain;
  domain.set_cubic(5);
  EXPECT_EQ(125, domain.volume());
}

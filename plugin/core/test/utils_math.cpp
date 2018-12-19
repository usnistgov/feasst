#include <gtest/gtest.h>
#include "core/include/utils_math.h"

TEST(UtilsMath, sgn) {
  EXPECT_EQ(feasst::sgn(-12.5), -1); 
  EXPECT_EQ(feasst::sgn(5./2.), 1); 
}

#include <gtest/gtest.h>
#include "core/include/utils_math.h"

namespace feasst {

TEST(UtilsMath, sgn) {
  EXPECT_EQ(sgn(-12.5), -1);
  EXPECT_EQ(sgn(5./2.), 1);
}

}  // namespace feasst

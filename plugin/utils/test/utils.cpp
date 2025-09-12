#include "utils/test/utils.h"
#include "utils/include/utils.h"

namespace feasst {
TEST(Utils, is_equal) {
  EXPECT_TRUE(is_equal(1.1, 1.1));
  std::vector<double> val = {1, 2, 3};
  EXPECT_TRUE(is_equal(val, val));
  std::vector<std::vector<double> > val2d;
  val2d.push_back({1, 2, 3});
  val2d.push_back({1, 2, 3});
  EXPECT_TRUE(is_equal(val2d, val2d));
  std::vector<std::vector<double> > vval2d;
  vval2d.push_back({1, 2, 3});
  vval2d.push_back({1, 2});
  EXPECT_FALSE(is_equal(val2d, vval2d));
}

}  // namespace feasst

#include <vector>
#include "utils/test/utils.h"
#include "utils/include/utils_io.h"
#include "flat_histogram/include/criteria_flat_histogram.h"

namespace feasst {

TEST(Macrostate, segment) {
  std::vector<double> windows = segment(0, 300, 12, 2);
  EXPECT_EQ(windows[3], 150);
  EXPECT_NEAR(windows[10], 273.861, 0.001);
}

TEST(Macrostate, windows) {
  std::vector<std::vector<int> > windows = window(0, 300, 12, 2);
  EXPECT_EQ(windows[3][0], 150);
  EXPECT_EQ(windows[3][1], 173);
}

}  // namespace feasst

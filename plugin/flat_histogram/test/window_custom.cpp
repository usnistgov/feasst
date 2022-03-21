#include <vector>
#include "utils/test/utils.h"
#include "flat_histogram/include/window_custom.h"

namespace feasst {

TEST(WindowCustom, window) {
  //WindowCustom windows({5, 100, 150, 200});
  WindowCustom windows({5, 100, 150, 200}, {{"overlap", "3"}});
  std::vector<double> segment = windows.segment();
  EXPECT_NEAR(segment[0], 5, 1e-13);
  EXPECT_NEAR(segment[1], 100, 1e-13);
  EXPECT_NEAR(segment[2], 150, 1e-13);
  EXPECT_NEAR(segment[3], 200, NEAR_ZERO);

  std::vector<std::vector<int> > bounds = windows.boundaries();
  EXPECT_EQ(bounds[0][0], 5);
  EXPECT_EQ(bounds[0][1], 100);
  EXPECT_EQ(bounds[1][0], 98);
  EXPECT_EQ(bounds[1][1], 150);
  EXPECT_EQ(bounds[2][0], 148);
  EXPECT_EQ(bounds[2][1], 200);
}

}  // namespace feasst

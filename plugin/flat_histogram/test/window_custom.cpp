#include <vector>
#include "utils/test/utils.h"
#include "flat_histogram/include/window_custom.h"

namespace feasst {

TEST(WindowCustom, window) {
  //WindowCustom windows({5, 100, 150, 200});
  WindowCustom win1({5, 98, 148, 200}, {{"overlap", "3"}});
  WindowCustom win2({{"min0", "5"}, {"min1", "98"}, {"min2", "148"}, {"max", "200"}, {"overlap", "3"}});
  for (const WindowCustom& win : {win1, win2}) {
    std::vector<double> segment = win.segment();
    EXPECT_EQ(segment.size(), 4);
    EXPECT_NEAR(segment[0], 5, 1e-13);
    EXPECT_NEAR(segment[1], 98, 1e-13);
    EXPECT_NEAR(segment[2], 148, 1e-13);
    EXPECT_NEAR(segment[3], 200, NEAR_ZERO);

    std::vector<std::vector<int> > bounds = win.boundaries();
    EXPECT_EQ(bounds.size(), 3);
    EXPECT_EQ(bounds[0][0], 5);
    EXPECT_EQ(bounds[0][1], 100);
    EXPECT_EQ(bounds[1][0], 98);
    EXPECT_EQ(bounds[1][1], 150);
    EXPECT_EQ(bounds[2][0], 148);
    EXPECT_EQ(bounds[2][1], 200);
  }
  WindowCustom win3({{"min0", "0"}, {"min1", "1"}, {"max", "1"}, {"overlap", "0"}});
  std::vector<std::vector<int> > bounds = win3.boundaries();
  EXPECT_EQ(bounds.size(), 2);
  EXPECT_EQ(bounds[0][0], 0);
  EXPECT_EQ(bounds[0][1], 0);
  EXPECT_EQ(bounds[1][0], 1);
  EXPECT_EQ(bounds[1][1], 1);
}

}  // namespace feasst

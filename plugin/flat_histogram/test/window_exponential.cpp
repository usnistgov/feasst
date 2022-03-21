#include <vector>
#include "utils/test/utils.h"
#include "flat_histogram/include/window_exponential.h"

namespace feasst {

TEST(WindowExponential, window) {
  WindowExponential windows({
    //{"minimum", "0"},
    {"maximum", "200"},
    {"num", "4"},
    {"alpha", "1.5"}});
  std::vector<double> segment = windows.segment();
  EXPECT_NEAR(segment[0], 0, NEAR_ZERO);
  EXPECT_NEAR(segment[1], 79.370052598409956, 1e-13);
  EXPECT_NEAR(segment[2], 125.99210498948727, 1e-13);
  EXPECT_NEAR(segment[3], 165.09636244473126, 1e-13);
  EXPECT_NEAR(segment[4], 200, NEAR_ZERO);

  std::vector<std::vector<int> > bounds = windows.boundaries();
  EXPECT_EQ(bounds[0][0], 0);
  EXPECT_EQ(bounds[0][1], 79);
  EXPECT_EQ(bounds[1][0], 79);
  EXPECT_EQ(bounds[1][1], 126);
  EXPECT_EQ(bounds[2][0], 126);
  EXPECT_EQ(bounds[2][1], 165);
  EXPECT_EQ(bounds[3][0], 165);
  EXPECT_EQ(bounds[3][1], 200);

  WindowExponential windows2({
    //{"minimum", "0"},
    {"maximum", "200"},
    {"num", "4"},
    {"alpha", "2"},
    {"overlap", "5"}});
  bounds = windows2.boundaries();
  EXPECT_EQ(bounds[0][0], 0);
  EXPECT_EQ(bounds[0][1], 100);
  EXPECT_EQ(bounds[1][0], 96);
  EXPECT_EQ(bounds[1][1], 141);
  EXPECT_EQ(bounds[2][0], 137);
  EXPECT_EQ(bounds[2][1], 173);
  EXPECT_EQ(bounds[3][0], 169);
  EXPECT_EQ(bounds[3][1], 200);

  TRY(
    WindowExponential({
      {"alpha", "1.75"},
      {"num", "32"},
      {"maximum", "56"},
      {"overlap", "3"}}).boundaries();
    CATCH_PHRASE("Windows are assumed to only overlap with the nearest");
  );
}

}  // namespace feasst

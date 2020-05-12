#include <vector>
#include "utils/test/utils.h"
//#include "utils/include/utils_io.h"
//#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/macrostate_num_particles.h"

namespace feasst {

// attempt to swap soft bounds
TEST(Macrostate, swap) {
  auto macro1 = MakeMacrostateNumParticles(
    Histogram({{"width", "1"}, {"max", "25"}}),
    {{"soft_max", "15"}, {"soft_min", "0"}}
  );
  auto macro2 = MakeMacrostateNumParticles(
    Histogram({{"width", "1"}, {"max", "25"}}),
    {{"soft_max", "25"}, {"soft_min", "10"}}
  );
  EXPECT_EQ(15, macro1->soft_max());
  EXPECT_EQ(0, macro1->soft_min());
  macro1->swap_soft_bounds(macro2.get());
  EXPECT_EQ(25, macro1->soft_max());
  EXPECT_EQ(10, macro1->soft_min());
}

}  // namespace feasst

#include <vector>
#include "utils/test/utils.h"
#include "flat_histogram/include/macrostate_energy.h"

namespace feasst {

// attempt to swap soft bounds
TEST(MacrostateEnergy, serialize) {
  auto macro = MakeMacrostateEnergy(Histogram({{"width", "1"}, {"max", "2"}}));
  test_serialize<MacrostateEnergy, Macrostate>(*macro);
}

}  // namespace feasst

#include <vector>
#include "utils/test/utils.h"
#include "flat_histogram/include/macrostate_position.h"

namespace feasst {

TEST(MacrostatePosition, serialize) {
  auto macro = MakeMacrostatePosition(Histogram({{"width", "1"}, {"max", "2"}}));
  test_serialize<MacrostatePosition, Macrostate>(*macro);
}

}  // namespace feasst

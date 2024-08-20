#include <vector>
#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "flat_histogram/include/macrostate_position.h"

namespace feasst {

TEST(MacrostatePosition, serialize) {
  auto macro = std::make_unique<MacrostatePosition>(Histogram({{"width", "1"}, {"max", "2"}}));
  test_serialize(macro);
}

}  // namespace feasst

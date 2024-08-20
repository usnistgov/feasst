#include <vector>
#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "model_expanded/include/macrostate_model.h"

namespace feasst {

TEST(MacrostateModel, serialize) {
  auto macro = std::make_unique<MacrostateModel>(Histogram({{"width", "1"}, {"max", "2"}}));
  test_serialize(macro);
}

}  // namespace feasst

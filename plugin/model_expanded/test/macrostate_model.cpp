#include <vector>
#include "utils/test/utils.h"
#include "model_expanded/include/macrostate_model.h"

namespace feasst {

TEST(MacrostateModel, serialize) {
  auto macro = MakeMacrostateModel(Histogram({{"width", "1"}, {"max", "2"}}));
  test_serialize<MacrostateModel, Macrostate>(*macro);
}

}  // namespace feasst

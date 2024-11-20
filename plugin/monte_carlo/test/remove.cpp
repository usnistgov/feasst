#include "utils/test/utils.h"
#include "monte_carlo/include/remove.h"

namespace feasst {

TEST(Remove, serialize) {
  auto remove = std::make_shared<Remove>(argtype({{"name", "something"}}));
  Action remove2 = test_serialize(*remove);
}

}  // namespace feasst

#include "utils/test/utils.h"
#include "actions/include/record.h"

namespace feasst {

TEST(Record, serialize) {
  auto obj = std::make_shared<Record>(argtype({{"save_positions", "something"}}));
  Action obj2 = test_serialize(*obj);
}

}  // namespace feasst

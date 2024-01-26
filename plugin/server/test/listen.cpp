#include "utils/test/utils.h"
#include "server/include/listen.h"

namespace feasst {

TEST(Listen, serialize) {
  auto obj = MakeListen();
  auto obj2 = test_serialize(*obj);
}

}  // namespace feasst

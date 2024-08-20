#include "utils/test/utils.h"
#include "server/include/server.h"
#include "server/include/listen.h"

namespace feasst {

TEST(Listen, serialize) {
  auto obj = std::make_unique<Listen>();
  auto obj2 = test_serialize(obj);
}

}  // namespace feasst

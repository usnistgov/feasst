#include "utils/test/utils.h"
#include "server/include/server.h"

namespace feasst {

TEST(Server, serialize) {
  auto obj = std::make_shared<Server>();
  auto obj2 = test_serialize(*obj);
}

}  // namespace feasst

#include "utils/test/utils.h"
#include "server/include/model_server.h"

namespace feasst {

TEST(ModelServer, serialize) {
  auto obj = std::make_unique<ModelServer>();
  auto obj2 = test_serialize(obj);
}

}  // namespace feasst

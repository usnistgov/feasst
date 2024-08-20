#include "utils/test/utils.h"
#include "server/include/visit_model_inner_server.h"
#include "server/include/server.h"

namespace feasst {

TEST(VisitModelInnerServer, serialize) {
  auto obj = std::make_unique<VisitModelInnerServer>();
  auto obj2 = test_serialize(obj);
}

}  // namespace feasst

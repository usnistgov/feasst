#include "utils/test/utils.h"
#include "server/include/visit_model_inner_server.h"

namespace feasst {

TEST(VisitModelInnerServer, serialize) {
  auto obj = MakeVisitModelInnerServer();
  auto obj2 = test_serialize(*obj);
}

}  // namespace feasst

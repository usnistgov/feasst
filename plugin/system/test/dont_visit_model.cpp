#include "utils/test/utils.h"
#include "system/include/dont_visit_model.h"

namespace feasst {

TEST(DontVisitModel, energy) {
  auto visit = MakeDontVisitModel();
  auto visit2 = test_serialize<DontVisitModel, VisitModel>(*visit);
}

}  // namespace feasst

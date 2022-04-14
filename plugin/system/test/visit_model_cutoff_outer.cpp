#include "utils/test/utils.h"
#include "system/include/visit_model_cutoff_outer.h"

namespace feasst {

TEST(VisitModelCutoffOuter, energy) {
  auto visit = MakeVisitModelCutoffOuter();
  auto visit2 = test_serialize<VisitModelCutoffOuter, VisitModel>(*visit);
}

}  // namespace feasst

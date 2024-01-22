#include <vector>
#include "utils/test/utils.h"
#include "model_expanded/include/perturb_model.h"

namespace feasst {

TEST(PerturbModel, serialize) {
  auto perturb = MakePerturbModel();
  PerturbModel perturb2 = test_serialize(*perturb);
}

}  // namespace feasst

#include "utils/test/utils.h"
#include "monte_carlo/include/perturb_add_remove.h"

namespace feasst {

TEST(PerturbAddRemove, serialize) {
  auto perturb = std::make_unique<PerturbAddRemove>();
  auto perturb2 = test_serialize_unique(*perturb);
}

}  // namespace feasst

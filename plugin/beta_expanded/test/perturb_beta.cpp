#include "utils/test/utils.h"
#include "beta_expanded/include/perturb_beta.h"

namespace feasst {

TEST(PerturbBeta, serialize) {
  auto beta = MakePerturbBeta({{"fixed_beta_change", "1"}});
  PerturbBeta beta2 = test_serialize(*beta);
}

}  // namespace feasst

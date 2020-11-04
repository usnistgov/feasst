#include "utils/test/utils.h"
#include "chain/include/perturb_branch.h"

namespace feasst {

TEST(PerturbBranch, serialize) {
  PerturbBranch branch;
  PerturbBranch branch2 = test_serialize(branch);
}

}  // namespace feasst

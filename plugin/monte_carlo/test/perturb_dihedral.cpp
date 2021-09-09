#include "utils/test/utils.h"
#include "monte_carlo/include/perturb_dihedral.h"

namespace feasst {

TEST(PerturbDihedral, serialize) {
  PerturbDihedral add;
  PerturbDihedral add2 = test_serialize(add);
}

}  // namespace feasst

#include "monte_carlo/test/perturb_test.h"

namespace feasst {

TEST(Perturb, Revert) {
  System system = default_system();
  test_revert(&system);
}

}  // namespace feasst

#include "utils/test/utils.h"
#include "monte_carlo/include/ref_potential.h"

namespace feasst {

TEST(RefPotential, serialize) {
  RefPotential obj;
  test_serialize(obj);
}

}  // namespace feasst

#include "utils/test/utils.h"
#include "monte_carlo/include/opt_potential.h"

namespace feasst {

TEST(OptPotential, serialize) {
  OptPotential obj;
  test_serialize(obj);
}

}  // namespace feasst

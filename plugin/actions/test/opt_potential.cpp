#include "utils/test/utils.h"
#include "actions/include/opt_potential.h"

namespace feasst {

TEST(OptPotential, serialize) {
  OptPotential obj;
  test_serialize(obj);
}

}  // namespace feasst

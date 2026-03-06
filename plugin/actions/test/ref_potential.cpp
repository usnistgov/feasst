#include "utils/test/utils.h"
#include "actions/include/ref_potential.h"

namespace feasst {

TEST(RefPotential, serialize) {
  RefPotential obj;
  test_serialize(obj);
}

}  // namespace feasst

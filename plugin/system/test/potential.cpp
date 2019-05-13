#include <sstream>
#include "utils/test/utils.h"
#include "system/include/potential.h"

namespace feasst {

TEST(Potential, serialize) {
  Potential potential;
  Potential potential2 = test_serialize(potential);
}

}  // namespace feasst

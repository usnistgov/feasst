#include "utils/test/utils.h"
#include "monte_carlo/include/perturb_distance_angle.h"

namespace feasst {

TEST(PerturbDistanceAngle, serialize) {
  PerturbDistanceAngle add;
  PerturbDistanceAngle add2 = test_serialize(add);
}

}  // namespace feasst

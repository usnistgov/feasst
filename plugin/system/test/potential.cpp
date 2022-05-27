#include <sstream>
#include "utils/test/utils.h"
#include "configuration/test/config_utils.h"
#include "system/include/potential.h"

namespace feasst {

TEST(Potential, serialize) {
  Potential potential;
  Potential potential2 = test_serialize(potential);
}

TEST(Potential, args) {
  Configuration config = lj_sample4();
  auto potential = MakePotential({{"Model", "LennardJones"}, {"sigma", "2"}});
  potential->precompute(&config);
  EXPECT_EQ(2, potential->model_params().select("sigma").value(0));
}

}  // namespace feasst

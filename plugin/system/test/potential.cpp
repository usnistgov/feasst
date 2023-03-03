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

TEST(Potential, model_params) {
  auto config = MakeConfiguration({{"cubic_box_length", "20"}, {"particle_type0", "../forcefield/spce.fstprt"}});
  auto potential = MakePotential({{"Model", "LennardJones"}, {"cutoff", "0"}, {"cutoff1", "5"}});
  potential->precompute(config.get());
  EXPECT_EQ(0, potential->model_params().select("cutoff").value(0));
  EXPECT_EQ(5, potential->model_params().select("cutoff").value(1));
}

}  // namespace feasst

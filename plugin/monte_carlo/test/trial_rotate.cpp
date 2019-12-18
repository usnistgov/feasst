#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/trial_rotate.h"
#include "configuration/test/configuration_test.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/metropolis.h"
#include "configuration/include/file_xyz.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

TEST(TrialRotate, spce) {
  System system;
  {
    auto config = MakeConfiguration(
      {{"particle_type", "../forcefield/data.spce"},
       {"cubic_box_length", "4"}
       });
    config->add_particle_of_type(0);
    system.add(*config);
  }

  system.add(Potential(MakeLennardJones()));

  Metropolis criteria;
  criteria.set_beta(1.0);
  auto rotate = MakeTrialRotate({{"tunable_param", "90."}, {"weight", "0.5"}});
  EXPECT_EQ(0.5, rotate->weight());
  FileXYZ file;
  file.write("tmp/before", system.configuration());
  RandomMT19937 random;
  rotate->attempt(&criteria, &system, &random);
  file.write("tmp/after", system.configuration());

  test_serialize(*rotate);
}

}  // namespace feasst

#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/domain.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"
#include "system/include/potential.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

TEST(TrialRotate, spce) {
  System system;
  system.add(MakeConfiguration({{"cubic_side_length", "20"},
    {"particle_type", "../particle/spce.fstprt"},
    {"add_particles_of_type0", "1"}}));
  system.add(MakePotential(MakeLennardJones()));
  system.set(MakeThermoParams({{"beta", "1.0"}}));
  Metropolis criteria;
  criteria.precompute(&system);
  criteria.set_current_energy(0);
  criteria.set_current_energy_profile({0.});
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

#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "configuration/include/file_xyz.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/analyze.h"
#include "chain/include/trial_grow_linear.h"

namespace feasst {

TEST(TrialGrowLinear, chain10) {
  System system;
  system.add(*MakeConfiguration({{"cubic_box_length", "12"},
    {"particle_type", "../forcefield/chain10.fstprt"},
    {"add_particles_of_type0", "1"}}));
  system.add(MakePotential(MakeLennardJones()));
  system.set(MakeThermoParams({{"beta", "100.0"}}));
  auto criteria = MakeMetropolis();
  auto grow = MakeTrialGrowLinear(
    std::make_shared<TrialComputeMove>(),
    {{"particle_type", "0"}}
  );
  TrialFactory factory;
  factory.add(grow);
  factory.precompute(criteria.get(), &system);
  FileXYZ file;
  file.write("tmp/before2", system.configuration());
  RandomMT19937 random;
  for (int i = 0; i < 50; ++i) {
    factory.attempt(criteria.get(), &system, &random);
    file.write("tmp/after2", system.configuration());
  }
  EXPECT_NE(0, system.configuration().particle(0).site(0).position().coord(0));

  test_serialize(*grow);
}

}  // namespace feasst

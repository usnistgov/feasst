#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "chain/include/trial_grow.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/metropolis.h"
#include "configuration/include/file_xyz.h"
#include "monte_carlo/include/analyze.h"
#include "chain/include/analyze_rigid_bonds.h"

namespace feasst {

TEST(TrialGrow, chain10) {
  System system;
  {
    Configuration config(MakeDomain({{"cubic_box_length", "12"}}),
                         {{"particle_type", "../forcefield/data.chain10"}});
    config.add_particle_of_type(0);
    system.add(config);
  }
  system.add(Potential(MakeLennardJones()));

  auto criteria = MakeMetropolis();
  criteria->set_beta(100.0);
  auto grow = MakeTrialGrowLinear(
    std::make_shared<TrialComputeMove>(),
    {{"particle_type", "0"}}
  );
  TrialFactory factory;
  factory.add(grow);
  factory.precompute(criteria.get(), &system);
  FileXYZ file;
  file.write("tmp/before2", system.configuration());
  AnalyzeRigidBonds checker;
  RandomMT19937 random;
  for (int i = 0; i < 50; ++i) {
    factory.attempt(criteria.get(), &system, &random);
    file.write("tmp/after2", system.configuration());
    checker.update(*criteria, system, factory);
  }
  EXPECT_NE(0, system.configuration().particle(0).site(0).position().coord(0));

  test_serialize(*grow);
}

}  // namespace feasst

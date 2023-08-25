#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "configuration/include/file_xyz.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/analyze.h"
#include "chain/include/trial_pivot.h"

namespace feasst {

TEST(TrialPivot, chain10) {
  System system;
  system.add(*MakeConfiguration({{"cubic_side_length", "12"},
    {"particle_type", "../particle/chain10.fstprt"},
    {"add_particles_of_type0", "1"}}));
  system.add(MakePotential(MakeLennardJones()));
  system.set(MakeThermoParams({{"beta", "100.0"}}));
  auto criteria = MakeMetropolis();
  auto pivot = MakeTrialPivot({{"tunable_param", "90."}});
  TrialFactory factory;
  factory.add(pivot);
  factory.precompute(criteria.get(), &system);
  FileXYZ file;
  file.write("tmp/before", system.configuration());
  RandomMT19937 random;
  for (int i = 0; i < 50; ++i) {
    factory.attempt(criteria.get(), &system, &random);
    file.write("tmp/after", system.configuration());
  }
  EXPECT_NE(0, system.configuration().particle(0).site(0).position().coord(0));

  test_serialize(*pivot);
}

}  // namespace feasst

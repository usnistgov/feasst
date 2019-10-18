#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "chain/include/trial_grow.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "system/include/system.h"
#include "system/include/model_lj.h"
#include "monte_carlo/include/criteria_metropolis.h"
#include "configuration/include/file_xyz.h"
#include "monte_carlo/include/analyze.h"
#include "chain/include/analyze_rigid_bonds.h"

namespace feasst {

TEST(TrialGrow, chain10) {
  System system;
  {
    Configuration config;
    config.set_domain(Domain().set_cubic(12));
    config.add_particle_type("../forcefield/data.chain10");
    config.add_particle_of_type(0);
    system.add(config);
  }
  { // add potentials to system
    Potential potential;
    potential.set_model(std::make_shared<ModelLJ>());
    potential.set_visit_model(std::make_shared<VisitModel>());
    system.add_to_unoptimized(potential);
  }

  auto criteria = std::make_shared<CriteriaMetropolis>();
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
    checker.update(criteria.get(), system, factory);
  }
  EXPECT_NE(0, system.configuration().particle(0).site(0).position().coord(0));

  test_serialize(*grow);
}

}  // namespace feasst

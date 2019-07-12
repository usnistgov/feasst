#include <gtest/gtest.h>
#include "chain/include/trial.h"
#include "monte_carlo/include/trial_factory.h"
#include "system/include/system.h"
#include "system/include/model_lj.h"
#include "monte_carlo/include/criteria_metropolis.h"
#include "configuration/include/file_xyz.h"
#include "monte_carlo/include/analyze.h"
#include "chain/include/analyze_rigid_bonds.h"

namespace feasst {

TEST(TrialPivot, chain10) {
  seed_random_by_date();
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
  auto pivot = MakeTrialPivot({{"tunable_param", "90."}});
  TrialFactory factory;
  factory.add(pivot);
  factory.precompute(criteria.get(), &system);
  FileXYZ file;
  file.write("tmp/before", system.configuration());
  AnalyzeRigidBonds checker;
  for (int i = 0; i < 50; ++i) {
    factory.attempt(criteria.get(), &system);
    file.write("tmp/after", system.configuration());
    checker.update(criteria.get(), system, factory);
  }
  EXPECT_NE(0, system.configuration().particle(0).site(0).position().coord(0));
}

}  // namespace feasst

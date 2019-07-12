#include <gtest/gtest.h>
#include "monte_carlo/include/trial.h"
#include "configuration/test/configuration_test.h"
#include "system/include/system.h"
#include "system/include/model_lj.h"
#include "monte_carlo/include/criteria_metropolis.h"
#include "configuration/include/file_xyz.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

TEST(TrialRotate, spce) {
  System system;
  {
    Configuration config;
    config.add_particle_type("../forcefield/data.spce");
    config.add_particle_of_type(0);
    system.add(config);
  }

  {
    Potential potential;
    potential.set_model(std::make_shared<ModelLJ>());
    potential.set_visit_model(std::make_shared<VisitModel>());
    system.add_to_unoptimized(potential);
  }

  CriteriaMetropolis criteria;
  criteria.set_beta(1.0);
  auto rotate = MakeTrialRotate({{"tunable_param", "90."}, {"weight", "0.5"}});
  EXPECT_EQ(0.5, rotate->weight());
  FileXYZ file;
  file.write("tmp/before", system.configuration());
  rotate->attempt(&criteria, &system);
  file.write("tmp/after", system.configuration());
}

}  // namespace feasst

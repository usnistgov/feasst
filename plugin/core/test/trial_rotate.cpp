#include <gtest/gtest.h>
#include "core/include/trial_translate.h"
#include "core/include/trial_factory.h"
#include "core/test/configuration_test.h"
#include "core/include/system.h"
#include "core/include/model_lj.h"
#include "core/include/criteria_metropolis.h"
#include "core/include/file_xyz.h"
#include "core/include/analyze.h"

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
  TrialRotate rotate;
  FileXYZ file;
  file.write("tmp/before", system.configuration());
  rotate.set_max_move(60);
  rotate.attempt(&criteria, &system);
  file.write("tmp/after", system.configuration());
}

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
  criteria->set_beta(1.0);
  auto pivot = std::make_shared<TrialPivot>();
  pivot->set_max_move(90);
  TrialFactory factory;
  factory.add(pivot);
  FileXYZ file;
  file.write("tmp/before", system.configuration());
  factory.attempt(criteria.get(), &system);
  file.write("tmp/after", system.configuration());
  RigidBondChecker checker;
  checker.update(criteria, system, factory);
}

}  // namespace feasst

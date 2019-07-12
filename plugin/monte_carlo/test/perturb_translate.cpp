#include <gtest/gtest.h>
#include "monte_carlo/include/perturb.h"

namespace feasst {

TEST(PerturbTranslate, position) {
  Configuration config1;
  config1.set_domain(Domain().set_cubic(5));
  config1.add_particle_type("../forcefield/data.atom");
  config1.add_particle_of_type(0);

  // add custom property to site to see if it reverts
  config1.add_site_property("banana", 1., 0, 0);

  System system;
  system.add(config1);
  // system.add_model(std::make_shared<ModelLJ>());
  const Configuration& config = system.configuration();
  EXPECT_NEAR(1., config.particle(0).site(0).property("banana"), NEAR_ZERO);
  for (int dim = 0; dim < config.dimension(); ++dim) {
    EXPECT_NEAR(0., config.particle(0).position().coord(dim), NEAR_ZERO);
    EXPECT_NEAR(0., config.particle(0).site(0).position().coord(dim), NEAR_ZERO);
  }
  PerturbTranslate perturb;
  TrialSelectParticleOfType tsel;
  tsel.select(&system);
  perturb.perturb(&system, &tsel);
  const int particle_index = tsel.mobile().particle_index(0);

  // did the particle actually move away from the origin?
  for (int dim = 0; dim < config.dimension(); ++dim) {
    EXPECT_NE(0., config.particle(particle_index).position().coord(dim));
    EXPECT_NE(0., config.particle(particle_index).site(0).position().coord(dim));
  }

  // change custom property (2.2) to see if revert puts it back to original (1).
  system.get_configuration()->set_site_property("banana", 2.2, 0, 0);

  EXPECT_NEAR(2.2, config.particle(0).site(0).property("banana"), NEAR_ZERO);
  perturb.revert(&system);

  // did the particle go back to the origin?
  for (int dim = 0; dim < config.dimension(); ++dim) {
    EXPECT_NEAR(0., config.particle(particle_index).position().coord(dim), NEAR_ZERO);
  }
  EXPECT_NEAR(1., config.particle(0).site(0).property("banana"), NEAR_ZERO);
}

}  // namespace feasst

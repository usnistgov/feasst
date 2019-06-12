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
  const std::vector<double> disp = {1.43, -2.5, 0.03};
  for (int dim = 0; dim < static_cast<int>(disp.size()); ++dim) {
    EXPECT_NEAR(0., config.particle(0).position().coord(dim), NEAR_ZERO);
  }
  Position trajectory(disp);
  PerturbTranslate perturb;
  TrialSelectParticleOfType tsel;
  tsel.select(&system);
  perturb.perturb(&system, &tsel);
  const int particle_index = tsel.mobile().particle_index(0);

  // change custom property to test revert
  system.get_configuration()->set_site_property("banana", 2.2, 0, 0);

  EXPECT_NEAR(2.2, config.particle(0).site(0).property("banana"), NEAR_ZERO);
  perturb.revert(&system);
  for (int dim = 0; dim < static_cast<int>(disp.size()); ++dim) {
    EXPECT_NEAR(0., config.particle(particle_index).position().coord(dim), NEAR_ZERO);
  }
  EXPECT_NEAR(1., config.particle(0).site(0).property("banana"), NEAR_ZERO);
}

}  // namespace feasst

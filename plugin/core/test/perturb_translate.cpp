#include <gtest/gtest.h>
#include "core/include/perturb_translate.h"

namespace feasst {

TEST(PerturbTranslate, position) {
  Configuration config1;
  config1.set_domain(Domain().set_cubic(5));
  config1.add_particle_type("../forcefield/data.atom");
  config1.add_particle(0);
  System system;
  system.add_configuration(config1);
  const Configuration& config = system.configuration();
  const std::vector<double> disp = {1.43, -2.5, 0.03};
  for (int dim = 0; dim < static_cast<int>(disp.size()); ++dim) {
    EXPECT_NEAR(0., config.particle(0).position().coord(dim), NEAR_ZERO);
  }
  Position trajectory(disp);
  PerturbTranslate perturb;
  perturb.select_random_particle(0, config);
  const int particle_index = perturb.selection().particle_index(0);
  INFO("pi " << particle_index);
  perturb.translate_selected_particle(trajectory, &system);
  for (int dim = 0; dim < static_cast<int>(disp.size()); ++dim) {
    EXPECT_NEAR(disp[dim], config.particle(particle_index).position().coord(dim), NEAR_ZERO);
  }
  perturb.revert();
  for (int dim = 0; dim < static_cast<int>(disp.size()); ++dim) {
    EXPECT_NEAR(0., config.particle(particle_index).position().coord(dim), NEAR_ZERO);
  }
}

}  // namespace feasst

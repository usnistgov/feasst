#include <gtest/gtest.h>
#include "core/include/perturb_translate.h"

TEST(PerturbTranslate, position) {
  feasst::Configuration config;
  config.set_domain(feasst::DomainCuboid().set_cubic(5));
  config.add_particle_type("../forcefield/data.atom");
  config.add_particle(0);
  feasst::System system;
  system.add_configuration(config);
  //system.default_system();
//  feasst::Configuration * config = system.configuration(0);
  const std::vector<double> disp = {1.43, -2.5, 0.03};
  for (int dim = 0; dim < static_cast<int>(disp.size()); ++dim) {
    EXPECT_NEAR(0., system.particle(0).position().coord(dim), 1e-15);
  }
  feasst::Position trajectory(disp);
  feasst::PerturbTranslate perturb;
  // const int particle_index = 0;
  // config->select_particle(particle_index);
  perturb.select_random_particle(0, config);
  const int particle_index = perturb.selection().particle_index(0);
  INFO("pi " << particle_index);
  perturb.translate_selected_particle(trajectory, &system);
  for (int dim = 0; dim < static_cast<int>(disp.size()); ++dim) {
    EXPECT_NEAR(disp[dim], system.particle(particle_index).position().coord(dim), 1e-15);
  }
  perturb.revert();
  for (int dim = 0; dim < static_cast<int>(disp.size()); ++dim) {
    EXPECT_NEAR(0., system.particle(particle_index).position().coord(dim), 1e-15);
  }
}

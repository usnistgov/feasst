#include <gtest/gtest.h>
#include "core/include/select.h"
#include "core/include/configuration.h"
#include "core/include/file_xyz.h"

TEST(Select, random_particle) {
  feasst::Configuration config;
  config.set_domain(feasst::DomainCuboid().set_cubic(7));
  config.add_particle_type("../forcefield/data.spce");
  config.add_particle(0);
  feasst::Select selected;
  selected.add_random_particle(config.particles());
  EXPECT_EQ(0, selected.particle_index(0));
  std::vector<int> sites = {0, 1, 2};
  EXPECT_EQ(sites, selected.site_indices(0));
}

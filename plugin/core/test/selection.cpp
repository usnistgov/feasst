#include <gtest/gtest.h>
#include "core/include/selection.h"
#include "core/include/configuration.h"
#include "core/include/file_xyz.h"

TEST(Selection, random_particle) {
  feasst::Configuration config;
  config.set_domain(feasst::DomainCuboid().set_cubic(7));
  config.add_particle_type("../forcefield/data.spce");
  config.add_particle(0);
  feasst::Selection selected;
  selected.add_random_particle(config.particles());
  EXPECT_EQ(0, selected.particle_index(0));
  std::vector<int> sites = {0, 1, 2};
  EXPECT_EQ(sites, selected.site_indices(0));
}

TEST(Selection, group) {
  feasst::Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  for (int part = 0; part < 100; ++part) {
    config.add_particle(0);
  }
  feasst::FileXYZ().load("../plugin/core/test/data/spce_sample_config_periodic1.xyz", &config);
  config.select(feasst::Group().add_site_type(0));
  EXPECT_EQ(100, config.selection().num_particles());
  EXPECT_EQ(100, config.selection().num_sites());
  EXPECT_EQ(0, config.selection().site_index(0, 0));
}

TEST(PositionSelection, load) {
  feasst::Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  for (int part = 0; part < 100; ++part) {
    config.add_particle(0);
  }
  feasst::FileXYZ().load("../plugin/core/test/data/spce_sample_config_periodic1.xyz", &config);
  config.select_last_particle();
  feasst::PositionSelection select(config.selection(), config.particles());
  EXPECT_EQ(1, select.particle_positions().size());
  EXPECT_EQ(3, select.site_positions()[0].size());
  EXPECT_NEAR(9.606656807740E+00, select.site_positions()[0][0].coord(2), 1e-15);
  EXPECT_NEAR(6.753621862080E+00, select.site_positions()[0][2].coord(0), 1e-15);
}

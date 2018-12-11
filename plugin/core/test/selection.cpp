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
  EXPECT_EQ(0, selected.selection()[0].first);
  std::vector<int> sites = {0, 1, 2};
  EXPECT_EQ(sites, selected.selection()[0].second);
//  for (const std::pair<int, std::vector<int> >& pair : selected.selection()) {
//    std::cout << pair.first << ": " << feasst::str(pair.second) << std::endl;
//  }
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
  EXPECT_EQ(0, config.selection().selection()[0].second[0]);
}


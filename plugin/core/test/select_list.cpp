#include <gtest/gtest.h>
#include "core/include/select_list.h"
#include "core/include/file_xyz.h"
#include "core/include/constants.h"

namespace feasst {

TEST(SelectList, group) {
  Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  for (int part = 0; part < 100; ++part) {
    config.add_particle_of_type(0);
  }
  FileXYZ().load("../plugin/core/test/data/spce_sample_config_periodic1.xyz", &config);
  SelectList select;

  // last
  select.last_particle_added(&config);
  EXPECT_EQ(1, select.particle_positions().size());
  EXPECT_EQ(3, select.site_positions()[0].size());
  EXPECT_NEAR(9.606656807740E+00, select.site_positions()[0][0].coord(2), NEAR_ZERO);
  EXPECT_NEAR(6.753621862080E+00, select.site_positions()[0][2].coord(0), NEAR_ZERO);

  // group
  select.clear();
  select.add(config, Group().add_site_type(0));
  EXPECT_EQ(100, select.num_particles());
  EXPECT_EQ(100, select.num_sites());
  EXPECT_EQ(0, select.site_index(0, 0));
}

}  // namespace feasst

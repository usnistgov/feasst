#include <gtest/gtest.h>
#include "system/include/select_list.h"
#include "configuration/include/file_xyz.h"
#include "math/include/constants.h"
#include "math/include/random_mt19937.h"

namespace feasst {

TEST(SelectList, group) {
  seed_random_by_date();
  Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  for (int part = 0; part < 100; ++part) {
    config.add_particle_of_type(0);
  }
  FileXYZ().load("../plugin/system/test/data/spce_sample_config_periodic1.xyz", &config);
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

  RandomMT19937 random;
  Select * selp = &select;
  Select sel2 = selp->random_particle(&random);
  // HWH DEBUG why cant the following work?
  // sel2 = select.random_particle();
  EXPECT_FALSE(select.is_equivalent(sel2));
  EXPECT_FALSE(select.str().empty());
  EXPECT_EQ(1, sel2.num_sites());
  EXPECT_EQ(0, sel2.site_indices(0)[0]);

  select.check();
}

//TEST(SelectList, random_end_segment_in_particle) {
//  seed_random_by_date();
//  Configuration config;
//  config.add_particle_type("../forcefield/data.chain10");
//  config.add_particle_of_type(0);
//  SelectList select;
//  select.random_end_segment_in_particle(0, config, 1);
//  EXPECT_EQ(1, select.num_sites());
//}

}  // namespace feasst

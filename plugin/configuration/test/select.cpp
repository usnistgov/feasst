#include <vector>
#include "utils/test/utils.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(Select, add_remove) {
  Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  config.add_particle_of_type(0);
  config.add_particle_of_type(0);
  Select oxygen;
  oxygen.add_site(0, 0);
  oxygen.add_site(1, 0);

  // individual sites of multiple particles
  EXPECT_EQ(oxygen.num_sites(), 2);
  EXPECT_EQ(oxygen.num_particles(), 2);
  std::vector<int> indices = {0, 1};
  EXPECT_EQ(oxygen.particle_indices(), indices);
  indices = {0};
  EXPECT_EQ(oxygen.site_indices(0), indices);
  EXPECT_EQ(oxygen.site_indices(1), indices);

  // individual particles
  Select part0;
  part0.add_particle(config.particle(0), 0);
  EXPECT_EQ(part0.num_sites(), 3);
  EXPECT_EQ(part0.num_particles(), 1);
  indices = {0};
  EXPECT_EQ(part0.particle_indices(), indices);
  indices = {0, 1, 2};
  EXPECT_EQ(part0.site_indices(0), indices);

  Select part1;
  part1.add_particle(config.particle(1), 1);
  EXPECT_EQ(part1.num_sites(), 3);
  EXPECT_EQ(part1.num_particles(), 1);
  indices = {1};
  EXPECT_EQ(part1.particle_indices(), indices);

  // add together
  Select all;
  all.add(part0);
  all.add(part1);
  EXPECT_EQ(all.num_sites(), 6);
  EXPECT_EQ(all.num_particles(), 2);
  indices = {0, 1};
  EXPECT_EQ(all.particle_indices(), indices);
  indices = {0, 1, 2};
  EXPECT_EQ(all.site_indices(0), indices);
  EXPECT_EQ(all.site_indices(1), indices);

  // adding already existing doesn't change select
  auto all_bak = all;
  all.add(part1);
  EXPECT_TRUE(all.is_equivalent(all_bak));

  // remove some sites of particles
  all.remove(oxygen);
  EXPECT_EQ(all.num_sites(), 4);
  EXPECT_EQ(all.num_particles(), 2);
  indices = {0, 1};
  EXPECT_EQ(all.particle_indices(), indices);
  indices = {1, 2};
  EXPECT_EQ(all.site_indices(0), indices);
  EXPECT_EQ(all.site_indices(1), indices);

  // serialize
  Select all2 = test_serialize(all);
  EXPECT_EQ(all2.site_indices(1), indices);
}

TEST(SelectGroup, serialize) {
  Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  config.add_particle_of_type(0);
  SelectGroup oxygen;
  oxygen.set_group(Group().add_site_type(1));
  SelectGroup oxygen2 = test_serialize(oxygen);
}

}  // namespace feasst

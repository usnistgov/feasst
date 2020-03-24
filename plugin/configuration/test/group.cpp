#include "utils/test/utils.h"
#include "configuration/include/group.h"
#include "configuration/include/file_lmp.h"

namespace feasst {

TEST(Group, remove_sites) {
  Particle particle = FileLMP().read("../forcefield/data.spce");
  EXPECT_EQ(3, particle.num_sites());

  Particle oxygen(particle);
  MakeGroup({{"add_site_type", "0"}})->remove_sites(&oxygen);
  EXPECT_EQ(1, oxygen.num_sites());

  feasst::Particle hydrogen(particle);
//  std::vector<int> full_to_partial, partial_to_full;
  MakeGroup({{"add_site_type", "1"}})->remove_sites(&hydrogen);
//                                                &full_to_partial,
//                                                &partial_to_full);
  EXPECT_EQ(2, hydrogen.num_sites());
//  EXPECT_EQ(3, static_cast<int>(full_to_partial.size()));
//  EXPECT_EQ(-1, full_to_partial[0]);
//  EXPECT_EQ(0, full_to_partial[1]);
//  EXPECT_EQ(1, full_to_partial[2]);
//  EXPECT_EQ(2, static_cast<int>(partial_to_full.size()));
//  EXPECT_EQ(1, partial_to_full[0]);
//  EXPECT_EQ(2, partial_to_full[1]);
}

TEST(Group, serialize) {
  auto grp = MakeGroup({{"add_site_type", "1"},
                          {"add_particle_type", "0"}});
  std::stringstream ss;
  grp->serialize(ss);
  Group grp2(ss);
  std::stringstream ss2;
  grp2.serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
  Particle particle = FileLMP().read("../forcefield/data.spce");
  EXPECT_EQ(grp->site_indices(particle), grp2.site_indices(particle));
  std::vector<int> indices = {1, 2};
  EXPECT_EQ(grp->site_indices(particle), indices);
}

}  // namespace feasst

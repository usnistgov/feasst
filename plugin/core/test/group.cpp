#include <gtest/gtest.h>
#include "core/include/group.h"
#include "core/include/file_lmp.h"

TEST(Group, remove_sites) {
  feasst::Particle particle = feasst::FileLMP().read("../forcefield/data.spce");
  EXPECT_EQ(3, particle.num_sites());

  feasst::Particle oxygen(particle);
  feasst::Group().add_site_type(0).remove_sites(&oxygen);
  EXPECT_EQ(1, oxygen.num_sites());

  feasst::Particle hydrogen(particle);
//  std::vector<int> full_to_partial, partial_to_full;
  feasst::Group().add_site_type(1).remove_sites(&hydrogen);
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

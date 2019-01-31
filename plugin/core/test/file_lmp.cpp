#include <gtest/gtest.h>
#include "core/include/file_lmp.h"
#include "core/include/file.h"
#include "core/include/debug.h"
#include "core/include/constants.h"

namespace feasst {

TEST(FileLMP, data_lj) {
  FileLMP lmp_file;
  Particle particle = lmp_file.read("../forcefield/data.lj");
  EXPECT_EQ(particle.num_sites(), 1);
  EXPECT_EQ(particle.site(0).type(), 0);
  Position position;
  for (double coord : particle.site(0).position().coord()) {
    EXPECT_NEAR(0., coord, NEAR_ZERO);
  }
  FileLMP().read_properties("../forcefield/data.lj", &particle);
  EXPECT_NEAR(1., particle.site(0).property("epsilon"), NEAR_ZERO);
  EXPECT_NEAR(1., particle.site(0).property("sigma"), NEAR_ZERO);
  EXPECT_NEAR(3., particle.site(0).property("cutoff"), NEAR_ZERO);
  try {
    particle.site(0).property("charge");
    CATCH_PHRASE("not found");
  }
  EXPECT_EQ(0, lmp_file.num_bonds());
  EXPECT_EQ(0, lmp_file.num_bond_types());
  EXPECT_EQ(0, lmp_file.num_angles());
  EXPECT_EQ(0, lmp_file.num_angle_types());
}

TEST(FileLMP, data_spce) {
  FileLMP lmp_file;
  Particle particle = lmp_file.read("../forcefield/data.spce");
  EXPECT_EQ(particle.num_sites(), 3);
  Site site;
  site = particle.site(0);
  EXPECT_EQ(site.type(), 0);
  EXPECT_NEAR(0., site.position().coord(0), NEAR_ZERO);
  EXPECT_NEAR(0., site.position().coord(1), NEAR_ZERO);
  EXPECT_NEAR(0., site.position().coord(2), NEAR_ZERO);
  site = particle.site(1);
  EXPECT_EQ(site.type(), 1);
  EXPECT_NEAR(1.0, site.position().coord(0), NEAR_ZERO);
  EXPECT_NEAR(0., site.position().coord(1), NEAR_ZERO);
  EXPECT_NEAR(0., site.position().coord(2), NEAR_ZERO);
  site = particle.site(2);
  EXPECT_EQ(site.type(), 1);
  EXPECT_NEAR(-0.333313247568237000, site.position().coord(0), NEAR_ZERO);
  EXPECT_NEAR(0.942816142731718000, site.position().coord(1), NEAR_ZERO);
  EXPECT_NEAR(0., site.position().coord(2), NEAR_ZERO);

  lmp_file.read_properties("../forcefield/data.spce", &particle);
  EXPECT_NEAR(0.650169581, particle.site(0).property("epsilon"), NEAR_ZERO);
  EXPECT_NEAR(3.16555789, particle.site(0).property("sigma"), NEAR_ZERO);
  EXPECT_NEAR(-0.8476, particle.site(0).property("charge"), NEAR_ZERO);
  EXPECT_NEAR(10., particle.site(0).property("cutoff"), NEAR_ZERO);
  EXPECT_NEAR(0, particle.site(1).property("epsilon"), NEAR_ZERO);
  EXPECT_NEAR(0, particle.site(1).property("sigma"), NEAR_ZERO);
  EXPECT_NEAR(0.4238, particle.site(1).property("charge"), NEAR_ZERO);
  EXPECT_NEAR(10., particle.site(1).property("cutoff"), NEAR_ZERO);

  // bonds
  DEBUG("bonds");
  EXPECT_EQ(2, lmp_file.num_bonds());
  EXPECT_EQ(1, lmp_file.num_bond_types());
  EXPECT_EQ(1, lmp_file.num_angles());
  EXPECT_EQ(1, lmp_file.num_angle_types());
  EXPECT_EQ(2, particle.num_bonds());
  try {
    particle.bond(0).property("doesnotexist");
    CATCH_PHRASE("not found");
  }
  EXPECT_NEAR(1., particle.bond(0).property("l0"), NEAR_ZERO);
  EXPECT_NEAR(450, particle.bond(0).property("k"), NEAR_ZERO);
}

}  // namespace feasst

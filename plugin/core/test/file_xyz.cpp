#include <gtest/gtest.h>
#include "core/include/file_xyz.h"

TEST(FileXYZ, load) {
  feasst::Configuration config;
  feasst::FileXYZ().load("../plugin/core/test/data/lj_sample_config_periodic4.xyz", &config);
  EXPECT_NEAR(config.domain().volume(), 512., 1e-15);
  EXPECT_EQ(30, config.num_particles());
  EXPECT_NEAR(config.particle(29).position().coord(1), 3.786335083587E+00, 1e-15);
  feasst::FileXYZ().write("tmp/print.xyz", config);
}

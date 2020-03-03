#include <gtest/gtest.h>
#include "configuration/include/file_xyz.h"
#include "math/include/constants.h"
#include "configuration/test/configuration_test.h"

namespace feasst {

TEST(FileXYZ, load) {
  Configuration config = lj_sample();
  config.check();
  EXPECT_NEAR(config.domain()->volume(), 512., NEAR_ZERO);
  EXPECT_EQ(30, config.num_particles());
  EXPECT_NEAR(config.particle(29).position().coord(1), 3.786335083587E+00, NEAR_ZERO);
  FileXYZ().write("tmp/print.xyz", config);
}

TEST(FileXYZ, load_spce) {
  Configuration config = spce_sample();
  config.check();
  EXPECT_NEAR(8000, config.domain()->volume(), NEAR_ZERO);
}

}  // namespace feasst

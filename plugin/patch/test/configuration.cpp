#include <gtest/gtest.h>
#include "core/test/configuration_test.h"
#include "core/include/configuration.h"

namespace feasst {

TEST(Configuration, patch) {
  Configuration config;
  config.add_particle_type("../plugin/patch/forcefield/data.patch_one");
  EXPECT_FALSE(config.unique_type(0).site(0).is_director());
  EXPECT_TRUE(config.unique_type(0).site(1).is_director());
}

}  // namespace feasst

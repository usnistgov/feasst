#include <gtest/gtest.h>
#include "core/include/perturb_configs.h"
#include "core/include/file_xyz.h"
#include "core/include/constants.h"

namespace feasst {

TEST(PerturbConfigs, transfer_particle) {
  System sys;
  Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  sys.add(config);
  sys.add(config);
  Configuration * donor_config = sys.get_configuration(0);
  Configuration * acceptor_config = sys.get_configuration(1);
  FileXYZ().load("../plugin/core/test/data/lj_sample_config_periodic4.xyz", donor_config);
  acceptor_config->set_side_length(donor_config->domain().side_length());
  EXPECT_EQ(30, donor_config->num_particles());
  EXPECT_EQ(0, acceptor_config->num_particles());
  EXPECT_EQ(8*8*8, donor_config->domain().volume());
  EXPECT_EQ(8*8*8, acceptor_config->domain().volume());
  const double x_prev = donor_config->particle(29).position().coord(2);

  PerturbConfigs perturb;
  for (int num = 0; num < 30; ++num) {
    perturb.transfer_particle(0, &sys, 0, 1);
    EXPECT_EQ(30 - 1 - num, donor_config->num_particles());
    EXPECT_EQ(num + 1, acceptor_config->num_particles());
  }
  EXPECT_NEAR(x_prev, acceptor_config->particle(29).position().coord(2), NEAR_ZERO);
}

}  // namespace feasst

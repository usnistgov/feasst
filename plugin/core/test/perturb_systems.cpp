#include <gtest/gtest.h>
#include "core/include/perturb_systems.h"
#include "core/include/file_xyz.h"
#include "core/include/constants.h"

namespace feasst {

TEST(PerturbSystems, transfer_particle) {
  System donor;
  donor.add_configuration(Configuration());
  Configuration * donor_config = donor.configuration(0);
  donor_config->add_particle_type("../forcefield/data.lj");
  auto acceptor = donor;
  Configuration * acceptor_config = acceptor.configuration(0);
  FileXYZ().load("../plugin/core/test/data/lj_sample_config_periodic4.xyz", donor_config);
  acceptor_config->set_domain(donor_config->domain());
  EXPECT_EQ(30, donor_config->num_particles());
  EXPECT_EQ(0, acceptor_config->num_particles());
  EXPECT_EQ(8*8*8, donor_config->domain().volume());
  EXPECT_EQ(8*8*8, acceptor_config->domain().volume());
  const double x_prev = donor_config->particle(29).position().coord(2);

  PerturbSystems perturb;
  for (int num = 0; num < 30; ++num) {
    perturb.transfer_particle(0, &donor, &acceptor);
    EXPECT_EQ(30 - 1 - num, donor_config->num_particles());
    EXPECT_EQ(num + 1, acceptor_config->num_particles());
  }
  EXPECT_NEAR(x_prev, acceptor_config->particle(29).position().coord(2), feasst::NEAR_ZERO);
}

}  // namespace feasst

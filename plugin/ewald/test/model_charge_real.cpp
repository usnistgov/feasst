#include <gtest/gtest.h>
#include "ewald/include/model_charge_real.h"
#include "core/include/file_xyz.h"

TEST(ModelChargeReal, SRSW_refconfig) {
  feasst::Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  EXPECT_EQ(pow(0.8476, 2), (*config.unique_types().model_params().mixed_charge())[0][0]);
  EXPECT_EQ(-0.8476*0.4238, (*config.unique_types().model_params().mixed_charge())[0][1]);
  EXPECT_EQ(-0.8476*0.4238, (*config.unique_types().model_params().mixed_charge())[1][0]);
  EXPECT_EQ(pow(0.4238, 2), (*config.unique_types().model_params().mixed_charge())[1][1]);
  for (int part = 0; part < 100; ++part) {
    config.add_particle(0);
  }
  feasst::FileXYZ().load("../plugin/core/test/data/spce_sample_config_periodic1.xyz", &config);
  // feasst::FileXYZ().print("tmp/printspce.xyz", config);
  feasst::ModelChargeReal model;
  model.set_alpha(5.6/config.domain().min_side_length());
  feasst::VisitModel visit;
  visit.loop_by_particle(config, model);
  EXPECT_NEAR(-4646.8607665992577, visit.energy(), 1e-11);
}

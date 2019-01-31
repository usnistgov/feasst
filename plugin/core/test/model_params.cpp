#include <string>
#include <gtest/gtest.h>
#include "core/include/configuration.h"

namespace feasst {

TEST(ModelParams, max) {
  feasst::Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  const ModelParams& model_params = config.unique_types().model_params();
  EXPECT_EQ(model_params.epsilon().mixed_max(), 0.650169581);
  EXPECT_EQ(model_params.charge().mixed_max(), 0.8476*0.8476);
  EXPECT_EQ(model_params.charge().max(), 0.4238);
  EXPECT_EQ(pow(0.8476, 2), config.unique_types().model_params().mixed_charge()[0][0]);
  EXPECT_EQ(-0.8476*0.4238, config.unique_types().model_params().mixed_charge()[0][1]);
  EXPECT_EQ(-0.8476*0.4238, config.unique_types().model_params().mixed_charge()[1][0]);
  EXPECT_EQ(pow(0.4238, 2), config.unique_types().model_params().mixed_charge()[1][1]);
  config.set_model_param("cutoff", 0, 5);
  EXPECT_EQ(model_params.cutoff().mixed_max(), 10);
  config.set_model_param("cutoff", 1, 5);
  EXPECT_EQ(model_params.cutoff().mixed_max(), 5);
}

}  // namespace feasst

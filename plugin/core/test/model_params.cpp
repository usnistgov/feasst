#include <string>
#include <gtest/gtest.h>
#include "core/include/configuration.h"

namespace feasst {

TEST(ModelParams, max) {
  feasst::Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  EXPECT_EQ(config.model_params().epsilon().mixed_max(), 0.650169581);
  EXPECT_EQ(config.model_params().charge().mixed_max(), 0.8476*0.8476);
  EXPECT_EQ(config.model_params().charge().max(), 0.4238);
  EXPECT_EQ(pow(0.8476, 2), config.model_params().mixed_charge()[0][0]);
  EXPECT_EQ(-0.8476*0.4238, config.model_params().mixed_charge()[0][1]);
  EXPECT_EQ(-0.8476*0.4238, config.model_params().mixed_charge()[1][0]);
  EXPECT_EQ(pow(0.4238, 2), config.model_params().mixed_charge()[1][1]);
  config.set_model_param("cutoff", 0, 5);
  EXPECT_EQ(config.model_params().cutoff().mixed_max(), 10);
  config.set_model_param("cutoff", 1, 5);
  EXPECT_EQ(config.model_params().cutoff().mixed_max(), 5);

  // custom model parameters
  ModelParams params = config.model_params();
  EXPECT_EQ(0.650169581, params.select("epsilon")->value(0));
  try {
    params.select("mistyped");
    CATCH_PHRASE("unrecognized name");
  }
  params.add(std::make_shared<ModelParam>());
  try {
    params.select("generic")->value(0);
    CATCH_PHRASE("size error");
  }
}

}  // namespace feasst

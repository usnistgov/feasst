#include <cmath>
#include <string>
#include "utils/test/utils.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(ModelParam, serialize) {
  ModelParam param;
  ModelParam param2 = test_serialize(param);
}

TEST(ModelParams, serialize) {
  ModelParams param;
  ModelParams param2 = test_serialize(param);
}

TEST(ModelParams, size) {
  feasst::Configuration config;
  config.add_particle_type("../particle/chain10.fstprt");
  EXPECT_EQ(1, config.model_params().size());
}

TEST(ModelParams, max) {
  feasst::Configuration config;
  config.add_particle_type("../particle/spce.fstprt");
  EXPECT_EQ(2, config.model_params().size());
  EXPECT_EQ(config.model_params().select("epsilon").mixed_max(), 0.650169581);
  EXPECT_EQ(config.model_params().select("charge").mixed_max(), 0.8476*0.8476);
  EXPECT_EQ(config.model_params().select("charge").max(), 0.4238);
  EXPECT_EQ(std::pow(0.8476, 2), config.model_params().select("charge").mixed_values()[0][0]);
  EXPECT_EQ(-0.8476*0.4238, config.model_params().select("charge").mixed_values()[0][1]);
  EXPECT_EQ(-0.8476*0.4238, config.model_params().select("charge").mixed_values()[1][0]);
  EXPECT_EQ(std::pow(0.4238, 2), config.model_params().select("charge").mixed_values()[1][1]);
  config.set_model_param("cutoff", 0, 5);
  EXPECT_EQ(config.model_params().select("cutoff").mixed_max(), 10);
  config.set_model_param("cutoff", 1, 5);
  EXPECT_EQ(config.model_params().select("cutoff").mixed_max(), 5);

  // custom model parameters
  ModelParams params = config.model_params();
  EXPECT_EQ(0.650169581, params.select("epsilon").value(0));
//  TRY(
//    params.select("mistyped");
//    CATCH_PHRASE("not found");
//  );
//  params.add(std::make_shared<ModelParam>());
//  TRY(
//    params.select("ModelParam").value(0);
//    CATCH_PHRASE("size error");
//  );

  // serialize
  ModelParams params2 = test_serialize(params);
  EXPECT_EQ(params.select("charge").mixed_values(), params2.select("charge").mixed_values());
  INFO(params2.str());
}

}  // namespace feasst

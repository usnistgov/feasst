#include "utils/test/utils.h"
#include "ewald/include/model_charge_self.h"
#include "configuration/test/configuration_test.h"

namespace feasst {

TEST(ModelChargeSelf, SRSW_refconfig) {
  Configuration config = spce_sample();
  config.add_model_param("alpha", 5.6/config.domain().min_side_length());
  ModelChargeSelf model;
  model.precompute(config.model_params());
  VisitModel visit;
  visit.compute(model, &config);
  EXPECT_NEAR(-23652.08040365018, visit.energy(), 1e-11);

  auto model2 = test_serialize<ModelChargeSelf, Model>(model);
  model2->compute(&config, &visit);
  EXPECT_NEAR(-23652.08040365018, visit.energy(), 1e-10);
}

}  // namespace feasst

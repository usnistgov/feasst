#include "utils/test/utils.h"
#include "ewald/include/charge_self.h"
#include "configuration/test/configuration_test.h"

namespace feasst {

TEST(ChargeSelf, SRSW_refconfig) {
  Configuration config = spce_sample();
  config.add_model_param("alpha", 5.6/config.domain().min_side_length());
  ChargeSelf model;
  model.precompute(config.model_params());
  VisitModel visit;
  visit.compute(model, &config);
  const double en = -23652.080370504391;
  EXPECT_NEAR(en, visit.energy(), 1e-11);

  auto model2 = test_serialize<ChargeSelf, Model>(model);
  model2->compute(&config, &visit);
  EXPECT_NEAR(en, visit.energy(), 1e-11);
}

}  // namespace feasst

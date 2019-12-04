#include "utils/test/utils.h"
#include "ewald/include/charge_screened.h"
#include "configuration/test/configuration_test.h"

namespace feasst {

TEST(ChargeScreened, SRSW_refconfig) {
  Configuration config = spce_sample();
  config.add_model_param("alpha", 5.6/config.domain().min_side_length());
  ChargeScreened model;
  model.precompute(config.model_params());
  VisitModel visit;
  visit.compute(model, &config);
  const double en = -4646.8607600872092;
  EXPECT_NEAR(en, visit.energy(), 1e-10);

  auto model2 = test_serialize<ChargeScreened, Model>(model);
  model2->compute(&config, &visit);
  EXPECT_NEAR(en, visit.energy(), 1e-10);
}

}  // namespace feasst

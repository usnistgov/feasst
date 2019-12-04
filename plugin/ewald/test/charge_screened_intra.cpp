#include "utils/test/utils.h"
#include "ewald/include/charge_screened_intra.h"
#include "configuration/test/configuration_test.h"
#include "system/include/visit_model_intra.h"

namespace feasst {

TEST(ChargeScreenedIntra, SRSW_refconfig) {
  Configuration config = spce_sample();
  config.add_model_param("alpha", 5.6/config.domain().min_side_length());
  ChargeScreenedIntra model;
  model.precompute(config.model_params());
  VisitModelIntra visit;
  visit.set_intra_cut(0);
  model.compute(&config, &visit);
  const double en = 23363.573741866534;
  EXPECT_NEAR(en, visit.energy(), 1e-11);

  auto model2 = test_serialize<ChargeScreenedIntra, Model>(model);
  model2->compute(&config, &visit);
  EXPECT_NEAR(en, visit.energy(), 1e-11);
}

}  // namespace feasst

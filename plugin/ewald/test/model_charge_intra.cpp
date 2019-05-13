#include "utils/test/utils.h"
#include "ewald/include/model_charge_intra.h"
#include "configuration/test/configuration_test.h"
#include "system/include/visit_model_intra.h"

namespace feasst {

TEST(ModelChargeIntra, SRSW_refconfig) {
  Configuration config = spce_sample();
  config.add_model_param("alpha", 5.6/config.domain().min_side_length());
  ModelChargeIntra model;
  model.precompute(config.model_params());
  VisitModelIntra visit;
  visit.set_intra_cut(0);
  model.compute(&config, &visit);
  EXPECT_NEAR(23363.573774608, visit.energy(), 1e-11);

  auto model2 = test_serialize<ModelChargeIntra, Model>(model);
  model2->compute(&config, &visit);
  EXPECT_NEAR(23363.573774608, visit.energy(), 1e-11);
}

}  // namespace feasst

#include <gtest/gtest.h>
#include "ewald/include/model_charge_intra.h"
#include "core/test/configuration_test.h"
#include "core/include/visit_model_intra.h"

namespace feasst {

TEST(ModelChargeIntra, SRSW_refconfig) {
  Configuration config = spce_sample();
  config.add_model_param("alpha", 5.6/config.domain().min_side_length());
  ModelChargeIntra model;
  VisitModelIntra visit;
  visit.set_intra_cut(0);
  visit.compute(model, &config, 0);
  EXPECT_NEAR(23363.573774608, visit.energy(), 1e-11);
}

}  // namespace feasst

#include <gtest/gtest.h>
#include "ewald/include/model_charge_screened.h"
#include "core/test/configuration_test.h"

namespace feasst {

TEST(ModelChargeScreened, SRSW_refconfig) {
  Configuration config = spce_sample();
  config.add_model_param("alpha", 5.6/config.domain().min_side_length());
  ModelChargeScreened model;
  VisitModel visit;
  visit.compute(model, &config);
  EXPECT_NEAR(-4646.8607665992467, visit.energy(), 1e-10);
}

}  // namespace feasst

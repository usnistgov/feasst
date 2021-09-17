#include "utils/test/utils.h"
#include "charge/include/charge_screened_intra.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_bond.h"
#include "charge/test/system_example.h"

namespace feasst {

TEST(ChargeScreenedIntra, SRSW_refconfig) {
  test_cases(
    { std::make_tuple(MakeCODATA2010(), 23363.573774608001),
      std::make_tuple(MakeCODATA2014(), 23363.573722131176),
      std::make_tuple(MakeCODATA2018(), 23363.573741866534) },
    MakeChargeScreenedIntra(),
    MakeVisitModelIntra({{"cutoff", "0"}})
  );

  test_cases(
    { std::make_tuple(MakeCODATA2010(), 23363.573774608001),
      std::make_tuple(MakeCODATA2014(), 23363.573722131176),
      std::make_tuple(MakeCODATA2018(), 23363.573741866534) },
    MakeChargeScreenedIntra(),
    MakeVisitModelBond()
  );
}

TEST(ChargeScreenedIntra, SRSW_refconfig_bond) {
  auto config = MakeConfiguration(MakeDomain({{"cubic_box_length", "20"}}),
                                  {{"particle_type", "../forcefield/spce.fstprt"}});
  config->add_model_param("alpha", 5.6/config->domain().min_side_length());
  config->add_particle_of_type(0);
  auto model = MakeChargeScreenedIntra();
  auto intra = MakeVisitModelIntra({{"cutoff", "0"}});
  auto bond = MakeVisitModelBond();
  Potential p1(model, intra);
  p1.precompute(config.get());
  const double p1en = p1.energy(config.get());
  Potential p2(model, bond);
  p2.precompute(config.get());
  EXPECT_NEAR(p1en, p2.energy(config.get()), NEAR_ZERO);
}

}  // namespace feasst

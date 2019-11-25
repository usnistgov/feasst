#include "utils/test/utils.h"
#include "system/include/visit_model.h"
#include "math/include/constants.h"
#include "system/test/system_test.h"

namespace feasst {

TEST(VisitModel, energy_map) {
  Configuration config = lj_sample();
  LennardJones model;
  VisitModel visit;
  visit.set_inner(MakeVisitModelInner(MakeEnergyMapAll()));
  visit.precompute(&config);
  model.compute(&config, &visit);
  const double en_lj_expect = -16.790321304625856;
  EXPECT_NEAR(en_lj_expect, visit.energy(), NEAR_ZERO);
  EXPECT_NEAR(en_lj_expect,
              visit.inner()->energy_map()->total(),
              10*NEAR_ZERO);
}

TEST(EnergyMapAll, revert) {
  auto map = MakeEnergyMapAll();
  map->revert();
}

}  // namespace feasst

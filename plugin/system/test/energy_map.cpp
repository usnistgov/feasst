#include "utils/test/utils.h"
#include "system/include/visit_model.h"
#include "math/include/constants.h"
#include "system/test/system_test.h"
#include "system/include/energy_map_all.h"

namespace feasst {

TEST(EnergyMap, energy_map) {
  Configuration config = lj_sample();
  LennardJones model;
  VisitModel visit;
  visit.set_inner(MakeVisitModelInner(MakeEnergyMapAll()));
  visit.precompute(&config);
  model.compute(&config, &visit);
  const double en_lj_expect = -16.790321304625856;
  EXPECT_NEAR(en_lj_expect, visit.energy(), NEAR_ZERO);
  EXPECT_NEAR(en_lj_expect,
              visit.inner()->energy_map()->total_energy(),
              1e-13);
}

}  // namespace feasst

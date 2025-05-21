#include <memory>
#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "system/test/sys_utils.h"

namespace feasst {

TEST(Criteria, current_energy) {
  System sys = two_particle_system();
  const double pe_expected = 4*(pow(1.25, -12) - pow(1.25, -6));
  EXPECT_NEAR(sys.energy(), pe_expected, NEAR_ZERO);
  Metropolis crit;
  auto trans = MakeTrialTranslate({{"tunable_param", "0.1"}});
  trans->precompute(&crit, &sys);
  crit.set_current_energy(sys.energy());
  crit.set_current_energy_profile({0.});
  crit.precompute(&sys);
  RandomMT19937 random;
  TRY(
    auto sys2 = sys;
    auto crit2 = crit;
    trans->attempt(&crit2, &sys2, &random);
    CATCH_PHRASE("must set ThermoParams");
  );
  sys.set(MakeThermoParams({{"beta", "1."}}));
  trans->attempt(&crit, &sys, &random);
  EXPECT_NEAR(sys.energy(), crit.current_energy(), NEAR_ZERO);
}

}  // namespace feasst

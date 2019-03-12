#include <memory>
#include <gtest/gtest.h>
#include "core/include/criteria_metropolis.h"
#include "core/include/trial_translate.h"
#include "core/test/system_test.h"

namespace feasst {

TEST(Criteria, running_energy) {
  System sys = default_system();
  const double pe_expected = 4*(pow(1.25, -12) - pow(1.25, -6));
  EXPECT_NEAR(sys.energy(), pe_expected, NEAR_ZERO);
  auto trans = MakeTrialTranslate({{"max_move", "0.1"}});
  CriteriaMetropolis crit;
  crit.set_running_energy(sys.energy());
  try {
    auto sys2 = sys;
    auto crit2 = crit;
    trans->attempt(&crit2, &sys2);
    CATCH_PHRASE("beta must be initialized");
  }
  crit.set_beta(1.);
  trans->attempt(&crit, &sys);
  EXPECT_NEAR(sys.energy(), crit.running_energy(), NEAR_ZERO);
}

}  // namespace feasst

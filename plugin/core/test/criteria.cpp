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
  TrialTranslate trans;
  CriteriaMetropolis crit;
  crit.set_running_energy(sys.energy());
  trans.attempt(&crit, &sys);
  EXPECT_NEAR(sys.energy(), crit.running_energy(), NEAR_ZERO);
}

}  // namespace feasst

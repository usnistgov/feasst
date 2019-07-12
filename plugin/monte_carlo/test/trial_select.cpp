#include <gtest/gtest.h>
#include "monte_carlo/include/trial_select.h"
#include "configuration/test/configuration_test.h"
#include "system/test/system_test.h"

namespace feasst {

TEST(TrialSelectParticle, select) {
  System system = default_system();
  EXPECT_EQ(2, system.configuration().num_particles());
  TrialSelectParticle select;
  TrialSelectParticleOfType select2;
  int num_zero = 0;
  int num_zero2 = 0;
  const int num_select = 1e2;
  for (int i = 0; i < num_select; ++i) {
    select.select(&system);
    select2.select(&system);
    EXPECT_EQ(1, select.mobile().num_sites());
    EXPECT_EQ(1, select2.mobile().num_sites());
    if (select.mobile().particle_index(0) == 0) ++num_zero;
    if (select2.mobile().particle_index(0) == 0) ++num_zero2;
  }
  EXPECT_NE(num_zero, num_select);
  EXPECT_NE(num_zero, 0);
  EXPECT_NE(num_zero2, num_select);
  EXPECT_NE(num_zero2, 0);
}

}  // namespace feasst

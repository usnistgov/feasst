#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "system/test/sys_utils.h"

namespace feasst {

TEST(TrialSelectParticle, select) {
  System system = two_particle_system();
  EXPECT_EQ(2, system.configuration().num_particles());
  TrialSelectParticle select;
  select.precompute(&system);
  int num_zero = 0;
  const int num_select = 1e2;
  RandomMT19937 random;
  for (int i = 0; i < num_select; ++i) {
    select.sel(&system, &random);
    EXPECT_EQ(1, select.mobile().num_sites());
    if (select.mobile().particle_index(0) == 0) ++num_zero;
  }
  EXPECT_NE(num_zero, num_select);
  EXPECT_NE(num_zero, 0);
}

}  // namespace feasst

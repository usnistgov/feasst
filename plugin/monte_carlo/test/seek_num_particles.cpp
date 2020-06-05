#include "utils/test/utils.h"
#include "system/include/utils.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

TEST(SeekNumParticles, seek) {
  MonteCarlo mc;
  mc.set(lennard_jones());
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  SeekNumParticles(50).with_trial_add().run(&mc);
  EXPECT_EQ(mc.configuration().num_particles(), 50);
  try {
    SeekNumParticles(50, {{"particle_type", "1"}}).with_trial_add().run(&mc);
    CATCH_PHRASE("type: 1 >= num");
  }
  try {
    SeekNumParticles(50).with_trial_add({{"particle_type", "1"}}).run(&mc);
    CATCH_PHRASE("add: 1 not equal to type of seek");
  }
}

}  // namespace feasst

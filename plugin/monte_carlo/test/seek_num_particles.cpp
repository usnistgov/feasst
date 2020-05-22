#include "utils/test/utils.h"
#include "monte_carlo/include/utils.h"
#include "monte_carlo/include/seek_num_particles.h"

namespace feasst {

TEST(SeekNumParticles, seek) {
  MonteCarlo mc;
  lennard_jones(&mc);
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

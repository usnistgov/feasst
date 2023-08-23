#include "utils/test/utils.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

TEST(SeekNumParticles, seek) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_side_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  SeekNumParticles(50).with_trial_add().run(&mc);
  EXPECT_EQ(mc.configuration().num_particles(), 50);
  TRY(
    SeekNumParticles(50, {{"particle_type", "1"}}).with_trial_add().run(&mc);
    CATCH_PHRASE("type: 1 >= num");
  );
  TRY(
    SeekNumParticles(50).with_trial_add({{"particle_type", "1"}}).run(&mc);
    CATCH_PHRASE("add: 1 not equal to type of seek");
  );
}

}  // namespace feasst

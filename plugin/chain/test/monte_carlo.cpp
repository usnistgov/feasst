#include <memory>
#include "utils/test/utils.h"
#include "utils/include/utils_io.h"
#include "math/include/accumulator.h"
#include "system/include/ideal_gas.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/check_physicality.h"
#include "steppers/include/check_properties.h"
#include "chain/include/analyze_rigid_bonds.h"
#include "chain/include/trial_grow.h"
#include "chain/include/trial_pivot.h"
#include "chain/include/trial_crankshaft.h"
#include "chain/include/trial_reptate.h"
#include "chain/include/trial_swap_sites.h"
#include "chain/include/recenter_particles.h"
#include "chain/test/system_chain.h"

namespace feasst {

TEST(MonteCarlo, chain) {
  MonteCarlo mc;
  mc.set(chain_system());
  mc.set(MakeMetropolis({{"beta", "1"}, {"chemical_potential", "1."}}));
  SeekNumParticles(1).with_trial_add().run(&mc);
  mc.add(MakeTrialTranslate({
    {"weight", "1."},
    {"tunable_param", "1."},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialRotate({
    {"weight", "1."},
    {"tunable_param", "20."},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialPivot({
    {"weight", "1."},
    {"tunable_param", "20."},
    {"max_length", "30"},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialReptate({
    {"weight", "1."},
    {"max_length", "1"},
//    {"reference_index", "0"},
//    {"num_steps", "2"},
  }));
  mc.add(MakeTrialCrankshaft({
    {"weight", "1."},
    {"tunable_param", "25."},
    {"max_length", "5."},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialGrowLinear(
    MakeTrialComputeMove(),
    {
//      {"weight", "0.1"},
      {"particle_type", "0"},
      {"num_steps", "3"},
      {"reference_index", "0"},
    }
  ));
  const int steps_per = 1e2;
  mc.add(MakeLog({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/chainlog.txt"},
  }));
  mc.add(MakeMovie({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/chain10movie.xyz"},
  }));
  mc.add(MakeCheckEnergy({
    {"steps_per", str(steps_per)},
    {"tolerance", "1e-10"},
  }));
  mc.add(MakeTuner({{"steps_per", str(steps_per)}}));
  mc.add(MakeAnalyzeRigidBonds({{"steps_per", str(steps_per)}}));
  mc.attempt(3e2);

  MonteCarlo mc2 = test_serialize(mc);
  EXPECT_EQ(mc2.analyzers().size(), 3);
}

}  // namespace feasst

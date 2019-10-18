#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "ewald/test/system_example.h"
#include "monte_carlo/include/criteria_metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"
#include "steppers/include/movie.h"
#include "steppers/include/log.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check_energy.h"

namespace feasst {

TEST(MonteCarlo, spce) {
  MonteCarlo mc;
  {
    System system;
    {
      Configuration config({
        {"cubic_box_length", "24.8586887"},
        {"particle_type", "../forcefield/data.spce"}
      });
      config.add_model_param("alpha", 5.6/config.domain().min_side_length());
      system.add(config);
      const int kmax_squared = 38;
      add_ewald_with(MakeModelLJ(), &config, &system, kmax_squared);
    }
    mc.set(system);
  }
  const int steps_per = 2e0;
  mc.set(MakeCriteriaMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeMovie({{"file_name", "tmp/spce.xyz"}}));
  mc.add(MakeLog({{"file_name", "tmp/spce_log.txt"}}));
  mc.add(MakeTuner({{"steps_per", str(steps_per)}}));
  mc.add(MakeCheckEnergy({{"tolerance", "1e-8"}, {"steps_per", str(steps_per)}}));

  // Theres something wrong with MC and Ewald
  // mc.seek_num_particles(1);
  // INFO(mc.system().configuration().num_particles());
  // mc.attempt(1e3);

  test_serialize(mc);
}

}  // namespace feasst

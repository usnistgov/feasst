#include <memory>
#include "utils/test/utils.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "monte_carlo/include/criteria_metropolis.h"
#include "utils/include/utils_io.h"
#include "math/include/accumulator.h"
#include "system/test/system_test.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"
#include "steppers/include/num_particles.h"
#include "pipeline/include/pipeline.h"

namespace feasst {

TEST(Pipeline, NVT_benchmark) {
  seed_random_by_date();
  seed_random();
  Pipeline mc;
  { System system;
    system.add(Configuration({
      {"cubic_box_length", "8"},
      {"particle_type0", "../forcefield/data.lj"}
    }));
    { Potential potential;
      potential.set_model(std::make_shared<ModelLJ>());
      potential.set_visit_model(std::make_shared<VisitModel>());
      system.add(potential);
    }
    mc.set(system);
  }
  crit_trial_analyze(&mc, 1);
  mc.set(MakeCriteriaMetropolis({{"beta", "1.2"}, {"chemical_potential", "-10."}}));
  mc.add(MakeTrialTranslate({{"weight", "10."}})); // add a second translate
  add_trial_transfer(&mc,{{"weight", "10."}, {"particle_type", "0"}});
//  mc.add(MakeTrialAdd({{"weight", "10."}, {"particle_type", "0"}}));
  // mc.seek_num_particles(50);
  // mc.attempt(1e6);  // ~3.5 seconds (now 4.1)
  // for now, pipeline isn't fully implemented
  //mc.attempt(1e4);
  // DEBUG("\n" << mc.timer_str());
}

}  // namespace feasst

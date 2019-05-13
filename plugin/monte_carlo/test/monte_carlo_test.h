#include <memory>
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/criteria_metropolis.h"
#include "utils/include/utils_io.h"
#include "math/include/accumulator.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_cell.h"
#include "system/include/model_lj.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check.h"

namespace feasst {

inline MonteCarlo mc_lj() {
  MonteCarlo mc;
//  const double cutoff = 2.;

  { System system;
    system.add(Configuration({
      {"cubic_box_length", "8"},
      {"particle_type0", "../forcefield/data.lj"},
    }));

    { Potential potential;
      potential.set_model(std::make_shared<ModelLJ>());
      potential.set_visit_model(std::make_shared<VisitModel>());
      potential.set_model_params(system.configuration());
//      potential.set_model_param("cutoff", 0, cutoff);
//      EXPECT_NEAR(potential.model_params().mixed_cutoff()[0][0], cutoff, NEAR_ZERO);
      system.add_to_unoptimized(potential); }

    { Potential lrc;
      lrc.set_visit_model(std::make_shared<LongRangeCorrections>());
      lrc.set_model_params(system.configuration());
//      lrc.set_model_param("cutoff", 0, cutoff);
//      EXPECT_NEAR(lrc.model_params().mixed_cutoff()[0][0], cutoff, NEAR_ZERO);
      system.add_to_unoptimized(lrc); }

    mc.set(system);
  }

  mc.set(MakeCriteriaMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"max_move", "1."}}));
  const int num_periodic = 1e4;
  mc.add(MakeLog(
   {{"steps_per", str(num_periodic)},
    {"file_name", "tmp/ljlog.txt"}}));
  mc.add(MakeMovie(
   {{"steps_per", str(num_periodic)},
    {"file_name", "tmp/ljmovie.xyz"}}));
  mc.add(MakeEnergyCheck(
   {{"steps_per", str(num_periodic)},
    {"tolerance", "1e-10"}}));
  mc.add(MakeTuner({{"steps_per", str(num_periodic)}}));
  return mc;
}

}  // namespace feasst

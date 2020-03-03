#include <memory>
#include "utils/include/utils_io.h"
#include "math/include/accumulator.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_cell.h"
#include "system/include/lennard_jones.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check.h"
#include "steppers/include/check_energy.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_select_bond.h"
#include "monte_carlo/include/trial_select_angle.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/perturb_distance_angle.h"

namespace feasst {

inline void crit_trial_analyze(MonteCarlo * mc,
    const int num_periodic = 1e4,
    const bool translate = true,
    const bool rotate = false) {
  mc->set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  if (translate) {
    mc->add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  }
  if (rotate) {
    mc->add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "1."}}));
  }
  mc->add(MakeLog(
   {{"steps_per", str(num_periodic)},
    {"file_name", "tmp/ljlog.txt"}}));
  mc->add(MakeMovie(
   {{"steps_per", str(num_periodic)},
    {"file_name", "tmp/ljmovie.xyz"}}));
  mc->add(MakeCheckEnergy(
   {{"steps_per", str(num_periodic)},
    {"tolerance", "1e-10"}}));
  mc->add(MakeTuner({{"steps_per", str(num_periodic)}}));
}

inline void mc_lj(MonteCarlo * mc,
    const double box_length = 8,
    const std::string data = "../forcefield/data.lj",
    const int num_periodic = 1e4,
    const bool translate = true,
    const bool lrc = true) {
//  const double cutoff = 2.;

  { System system;
    system.add(Configuration(
      MakeDomain({{"cubic_box_length", str(box_length)}}),
      {{"particle_type0", data},
      //{"init_cells", "1."}
    }));
    system.add(Potential(MakeLennardJones()));
//    { Potential potential(MakeLennardJones());
//      potential.set_model_params(system.configuration());
////      potential.set_model_param("cutoff", 0, cutoff);
////      EXPECT_NEAR(potential.model_params().mixed_cutoff()[0][0], cutoff, NEAR_ZERO);
//      system.add_to_unoptimized(potential); }
    if (lrc) {
      system.add(Potential(MakeLongRangeCorrections()));
    }
//    { Potential lrc(MakeLongRangeCorrections());
//      lrc.set_model_params(system.configuration());
////      lrc.set_model_param("cutoff", 0, cutoff);
////      EXPECT_NEAR(lrc.model_params().mixed_cutoff()[0][0], cutoff, NEAR_ZERO);
//      //system.add_to_unoptimized(lrc);
//    }

    mc->set(system);
  }
  crit_trial_analyze(mc, num_periodic, translate);
}

// flag == 0, move; 1, add; 2, remove
inline std::shared_ptr<Trial> build_(const int flag, const std::string data) {
  auto trial = MakeTrial({{"weight", "100"}});
  auto select = MakeTrialSelectParticle({
    {"particle_type", "0"},
    {"site", "0"},
  });
  std::shared_ptr<Perturb> perturb = std::make_shared<PerturbAnywhere>();
  if (flag == 1) perturb = std::make_shared<PerturbAdd>();
  if (flag == 2) perturb = std::make_shared<PerturbRemove>();
  trial->add_stage(
    select,
    perturb,
    {{"num_steps", "4"}}
  );
  trial->add_stage(
    MakeTrialSelectBond({
      {"particle_type", "0"},
      {"mobile_site", "1"},
      {"anchor_site", "0"}
    }),
    std::make_shared<PerturbDistance>(),
    {{"num_steps", "4"}}
  );
  if (data == "../forcefield/data.spce") {
    trial->add_stage(
      MakeTrialSelectAngle({
        {"particle_type", "0"},
        {"mobile_site", "2"},
        {"anchor_site", "0"},
        {"anchor_site2", "1"}
      }),
      std::make_shared<PerturbDistanceAngle>(),
      {{"num_steps", "4"}}
    );
  }
  if (flag == 0) {
    trial->set(std::make_shared<TrialComputeMove>());
  } else if (flag == 1) {
    trial->set(std::make_shared<TrialComputeAdd>());
  } else if (flag == 2) {
    trial->set(std::make_shared<TrialComputeRemove>());
  } else {
    ERROR("unrecognized flag: " << flag);
  }
  return trial;
}

}  // namespace feasst

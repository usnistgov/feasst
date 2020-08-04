#include <memory>
#include "utils/include/io.h"
#include "math/include/accumulator.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_cell.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/trial_select_bond.h"
#include "monte_carlo/include/trial_select_angle.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/perturb_distance_angle.h"

namespace feasst {

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
  if (data == "forcefield/data.spce") {
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

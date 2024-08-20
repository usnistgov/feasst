#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialAnywhere(argtype args) {
  auto trial = MakeTrialMove(std::make_shared<TrialSelectParticle>(&args),
    std::make_shared<PerturbAnywhere>(),
    "TrialAnywhere",
    &args);
  feasst_check_all_used(args);
  return trial;
}

}  // namespace feasst

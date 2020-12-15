#include "utils/include/debug.h"
#include "monte_carlo/include/perturb_move.h"

namespace feasst {

void PerturbMove::perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held) {
  if (is_position_held) {
    select->set_trial_state(0);
    return;
  }
  move(system, select, random);
  set_revert_possible(true, select);
  set_finalize_possible(true, select);
  select->set_trial_state(1);
}

void PerturbMove::revert(System * system) {
  DEBUG("revert? " << revert_possible());
  if (revert_possible()) {
    Configuration* config = system->get_configuration();
    config->update_positions(revert_select()->mobile_original(),
      // don't wrap if reverting
      false);
    DEBUG("mobile orig " << revert_select()->mobile_original().str());
    system->revert(revert_select()->mobile_original());
  }
}

void PerturbMove::finalize(System * system) {
  if (finalize_possible()) {
    system->finalize(finalize_select()->mobile());
  }
}

}  // namespace feasst

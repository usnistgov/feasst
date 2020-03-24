#include "utils/include/debug.h"
#include "monte_carlo/include/acceptance.h"

namespace feasst {

void Acceptance::reset() {
  set_ln_metropolis_prob();
  set_reject();
  energy_new_ = 0.;
  energy_old_ = 0.;
  macrostate_shift_ = 0;
  macrostate_shift_type_ = 0;
  perturbed_.clear();
}

void Acceptance::add_to_perturbed(const Select& select) {
  perturbed_.add(select);
}

void Acceptance::set_perturbed_state(const int state) {
  DEBUG("state " << state);
  perturbed_.set_trial_state(state);
}

}  // namespace feasst


#include "monte_carlo/include/acceptance.h"

namespace feasst {

void Acceptance::reset() {
  set_ln_metropolis_prob();
  set_reject();
  energy_new_ = 0.;
  energy_old_ = 0.;
  macrostate_shift_ = 0;
  perturbed_.clear();
}

}  // namespace feasst

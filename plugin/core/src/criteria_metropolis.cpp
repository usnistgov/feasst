
#include "core/include/criteria_metropolis.h"

namespace feasst {

bool CriteriaMetropolis::is_accepted(const AcceptanceCriteria accept_criteria) {
  if ( (accept_criteria.force_rejection != 1) &&
       (random_.uniform() < exp(accept_criteria.ln_metropolis_prob)) ) {
    set_running_energy(accept_criteria.energy_new);
    return true;
  }
  return false;
}

}  // namespace feasst

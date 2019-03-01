
#include "flat_histogram/include/criteria_flat_histogram.h"
#include "core/include/utils_io.h"

namespace feasst {

bool CriteriaFlatHistogram::is_accepted(
    const AcceptanceCriteria accept_criteria) {
  ASSERT(bias_ != NULL, "bias must be initialized before trials");
  bool is_accepted;
  double ln_metropolis_prob = accept_criteria.ln_metropolis_prob;
  if (verbose) {
    cout << "accepted? " << exp(accept_criteria.ln_metropolis_prob +
                                bias_->ln_bias(macrostate_new_,
                                               macrostate_old_)) << endl;
  }
  if (accept_criteria.force_rejection == 1 ||
      !macrostate_->is_in_range(accept_criteria.system, this)) {
    is_accepted = false;
    ln_metropolis_prob = 0.;
    if (verbose) cout << "forced rejection" << endl;
  } else {
    after_attempt_(accept_criteria.system);
    if (random_.uniform() < exp(accept_criteria.ln_metropolis_prob +
                                       bias_->ln_bias(macrostate_new_,
                                                     macrostate_old_))) {
      is_accepted = true;
      if (verbose) cout << "accept" << endl;
    } else {
      is_accepted = false;
      if (verbose) cout << "reject" << endl;
    }
    if (verbose) {
      cout << "macro old new " << macrostate_old_ << " "
           << macrostate_new_ << endl;
    }
  }
  bias_->update(macrostate_old_,
                macrostate_new_,
                ln_metropolis_prob,
                is_accepted);
  if (is_accepted) {
    set_running_energy(accept_criteria.energy_new);
  }
  return is_accepted;
}

}  // namespace feasst

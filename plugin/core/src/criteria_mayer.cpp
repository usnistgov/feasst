
#include "core/include/criteria_mayer.h"
#include "core/include/utils_io.h"

namespace feasst {

bool CriteriaMayer::is_accepted(const AcceptanceCriteria accept_criteria) {
  const double energy_new = accept_criteria.energy_new;
  const double kF12 = exp(-beta()*energy_new) - 1.;
  bool is_accepted;
  if (verbose) cout << "energy new " << energy_new << " f12 " << kF12 << endl;
  if ( (accept_criteria.force_rejection != 1) &&
       (random_.uniform() < std::abs(kF12)/std::abs(f12old_)) ) {
    set_running_energy(energy_new);
    f12old_ = kF12;
    is_accepted = true;
    if (verbose) cout << "computing ref" << endl;
    const double energy_reference =
      accept_criteria.system->reference_energy(reference_index_);
    f12ref_ = exp(-beta()*energy_reference) - 1.;
    if (verbose) cout << "f12ref " << f12ref_ << endl;
  } else {
    is_accepted = false;
  }
  if (f12old_ < 0) {
    mayer_.accumulate(-1.);
  } else {
    mayer_.accumulate(1.);
  }
  mayer_ref_.accumulate(f12ref_/std::abs(f12old_));
  if (verbose) cout << "is accepted? " << is_accepted << endl;
  return is_accepted;
}

}  // namespace feasst

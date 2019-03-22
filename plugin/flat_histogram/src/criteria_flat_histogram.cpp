
#include "flat_histogram/include/criteria_flat_histogram.h"
#include "core/include/utils_io.h"

namespace feasst {

bool CriteriaFlatHistogram::is_accepted(
    const AcceptanceCriteria accept_criteria) {
  ASSERT(bias_ != NULL, "bias must be initialized before trials");
  bool is_accepted;
  double ln_metropolis_prob = accept_criteria.ln_metropolis_prob;
//  DEBUG("accepted? " << exp(accept_criteria.ln_metropolis_prob +
//                        bias_->ln_bias(macrostate_new_,
//                                       macrostate_old_)));
  if (accept_criteria.force_rejection == 1 ||
      !macrostate_->is_in_range(accept_criteria.system, this)) {
    is_accepted = false;
    ln_metropolis_prob = -NEAR_INFINITY;
    macrostate_new_ = macrostate_old_;
    DEBUG("forced rejection");
  } else {
    after_attempt_(accept_criteria.system);
    DEBUG("bias " << bias_->ln_bias(macrostate_new_, macrostate_old_)
      << " old " << macrostate_old_ << " new " << macrostate_new_);
    DEBUG("ln new " << bias_->ln_macro_prob().value(macrostate_new_));
    DEBUG("ln old " << bias_->ln_macro_prob().value(macrostate_old_));
    DEBUG("ln met " << ln_metropolis_prob);
    DEBUG("ln tot " << ln_metropolis_prob + bias_->ln_bias(macrostate_new_, macrostate_old_));
    if (random_.uniform() < exp(ln_metropolis_prob +
                                bias_->ln_bias(macrostate_new_,
                                               macrostate_old_))) {
      is_accepted = true;
      DEBUG("accept");
    } else {
      is_accepted = false;
      DEBUG("reject");
    }
    DEBUG("macro old new " << macrostate_old_ << " " << macrostate_new_);
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

std::string CriteriaFlatHistogram::write() const {
  std::stringstream ss;
  ss << Criteria::write();
  ss << bias_->write();
  ss << "macrostate "
     << bias_->write_per_bin_header() << " "
     << bin_trackers_.write_per_bin_header() << std::endl;
  const Histogram& hist = macrostate_->histogram();
  for (int bin = 0; bin < hist.size(); ++bin) {
    ss << hist.center_of_bin(bin) << " "
       << bias_->write_per_bin(bin) << " "
       << bin_trackers_.write_per_bin(bin) << std::endl;
  }
  return ss.str();
}

}  // namespace feasst

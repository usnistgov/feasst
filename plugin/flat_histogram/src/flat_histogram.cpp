
#include "flat_histogram/include/flat_histogram.h"
#include "utils/include/utils_io.h"

namespace feasst {

bool FlatHistogram::is_accepted(const Acceptance& acceptance,
    const System * system,
    const double uniform_random) {
  ASSERT(bias_ != NULL, "bias must be initialized before trials");
  bool is_accepted;
  double ln_metropolis_prob = acceptance.ln_metropolis_prob();
  DEBUG("macroshift " << acceptance.macrostate_shift());
  const int shift = acceptance.macrostate_shift()*num_trial_states();
  if (acceptance.reject() or
      !macrostate_->is_allowed(system, this, shift)) {
    is_accepted = false;
    ln_metropolis_prob = -NEAR_INFINITY;
    macrostate_new_ = macrostate_old_;
    DEBUG("forced rejection");
  } else {
    // the shift factory multiplied by number of states assumes only one
    // particle is added during an entire growth expanded cycle
    macrostate_new_ = macrostate_->bin(system, this) + shift;
    DEBUG("old " << macrostate_old_ << " new " << macrostate_new_);
    DEBUG("bias " << bias_->ln_bias(macrostate_new_, macrostate_old_));
    DEBUG("ln new " << bias_->ln_prob().value(macrostate_new_));
    DEBUG("ln old " << bias_->ln_prob().value(macrostate_old_));
    DEBUG("ln met " << ln_metropolis_prob);
    DEBUG("ln tot " << ln_metropolis_prob + bias_->ln_bias(macrostate_new_, macrostate_old_));
    if (uniform_random < exp(ln_metropolis_prob +
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
    set_current_energy(acceptance.energy_new());
    DEBUG("current energy: " << current_energy());
  } else {
    // return the macrostate to the current value, as used by Analyze, etc.
    macrostate_new_ = macrostate_old_;
  }
  was_accepted_ = is_accepted;
  last_acceptance_ = acceptance;
  return is_accepted;
}

std::string FlatHistogram::write() const {
  std::stringstream ss;
  ss << Criteria::write();
  ss << bias_->write();
  ss << "macrostate "
     << bias_->write_per_bin_header() << " "
     << std::endl;
  const Histogram& hist = macrostate_->histogram();
  for (int bin = 0; bin < hist.size(); ++bin) {
    ss << hist.center_of_bin(bin) << " "
       << bias_->write_per_bin(bin) << " "
       << std::endl;
  }
  return ss.str();
}

void FlatHistogram::revert(const bool accepted, const double ln_prob) {
  Criteria::revert(accepted, ln_prob);
  bias_->update_or_revert(macrostate_old_,
                          macrostate_new_,
                          ln_prob,
                          accepted,
                          true);
  macrostate_new_ = macrostate_old_;
}

class MapFlatHistogram {
 public:
  MapFlatHistogram() {
    FlatHistogram().deserialize_map()["FlatHistogram"] =
      MakeFlatHistogram();
  }
};

static MapFlatHistogram mapper_ = MapFlatHistogram();

FlatHistogram::FlatHistogram(std::istream& istr)
  : Criteria(istr) {
  feasst_deserialize_version(istr);
  // feasst_deserialize_fstdr(bias_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      bias_ = bias_->deserialize(istr);
    }
  }
  // feasst_deserialize_fstdr(macrostate_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      macrostate_ = macrostate_->deserialize(istr);
    }
  }
  feasst_deserialize(&macrostate_old_, istr);
  feasst_deserialize(&macrostate_new_, istr);
  feasst_deserialize(&is_macrostate_set_, istr);
}

void FlatHistogram::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(1, ostr);
  feasst_serialize_fstdr(bias_, ostr);
  feasst_serialize_fstdr(macrostate_, ostr);
  feasst_serialize(macrostate_old_, ostr);
  feasst_serialize(macrostate_new_, ostr);
  feasst_serialize(is_macrostate_set_, ostr);
}

}  // namespace feasst

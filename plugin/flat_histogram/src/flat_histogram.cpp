#include <cmath>
#include "utils/include/serialize.h"
#include "flat_histogram/include/flat_histogram.h"
#include "math/include/constants.h"

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
      !macrostate_->is_allowed(system, this, acceptance)) {
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
    macrostate_current_ = macrostate_new_;
  } else {
    // return the macrostate to the current value, as used by Analyze, etc.
    macrostate_current_ = macrostate_old_;
  }
  was_accepted_ = is_accepted;
  last_acceptance_ = acceptance;
  return is_accepted;
}

std::string FlatHistogram::write() const {
  std::stringstream ss;
  ss << Criteria::write();
  ss << bias_->write();
  ss << "macrostate,"
     << bias_->write_per_bin_header()
     << std::endl;
  const Histogram& hist = macrostate_->histogram();
  for (int bin = 0; bin < hist.size(); ++bin) {
    ss << hist.center_of_bin(bin) << ","
       << bias_->write_per_bin(bin)
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
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 937, "version mismatch: " << version);
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
  feasst_deserialize(&macrostate_current_, istr);
  feasst_deserialize(&is_macrostate_set_, istr);
}

void FlatHistogram::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(937, ostr);
  feasst_serialize_fstdr(bias_, ostr);
  feasst_serialize_fstdr(macrostate_, ostr);
  feasst_serialize(macrostate_old_, ostr);
  feasst_serialize(macrostate_new_, ostr);
  feasst_serialize(macrostate_current_, ostr);
  feasst_serialize(is_macrostate_set_, ostr);
}

double FlatHistogram::pressure(const double volume, const int phase) const {
  int min, max;
  phase_boundary_(phase, &min, &max);
  return (-ln_prob_().value(0) + std::log(ln_prob_().sum_probability(min, max)))
         /volume/beta();
}

double FlatHistogram::average(const LnProbability& ln_prob,
     const std::vector<double>& macrostate_averages,
     const int phase) const {
  ASSERT(ln_prob.size() == static_cast<int>(macrostate_averages.size()),
    "size mismatch: ln_prob:" << ln_prob.size() <<
    " macro:" << macrostate_averages.size());
  int min, max;
  phase_boundary_(ln_prob, phase, &min, &max);
  double average = 0.;
  for (int bin = min; bin < max + 1; ++bin) {
    average += macrostate_averages[bin]*std::exp(ln_prob.value(bin));
  }
  return average/ln_prob.sum_probability(min, max);
}

double FlatHistogram::average_macrostate(const LnProbability& ln_prob,
    const int phase) const {
  int min, max;
  phase_boundary_(phase, &min, &max);
  double average = 0.;
  for (int bin = min; bin < max + 1; ++bin) {
    average += macrostate_->value(bin)*std::exp(ln_prob.value(bin));
  }
  return average/ln_prob.sum_probability(min, max);
}

void FlatHistogram::set(const std::shared_ptr<Bias> bias) {
  ASSERT(is_macrostate_set_, "set macrostate before bias");
  bias_ = bias;
  bias_->resize(macrostate_->histogram());
}

void FlatHistogram::before_attempt(const System* system) {
  macrostate_old_ = macrostate_->bin(system, this);
  DEBUG("macro old " << macrostate_old_);
}

LnProbability FlatHistogram::reweight(const double delta_conjugate) {
  LnProbability lnpirw = deep_copy(bias()->ln_prob());
  for (int macro = 0; macro < lnpirw.size(); ++macro) {
    lnpirw.add(macro, macrostate()->histogram().center_of_bin(macro)
               *delta_conjugate);
  }
  lnpirw.normalize();
  return lnpirw;
}

void FlatHistogram::phase_boundary_(const LnProbability& ln_prob,
    const int phase, int * min, int * max) const {
  std::vector<int> mins = ln_prob.minima();
  const double num_min = static_cast<int>(mins.size());
  if (num_min == 0) {
    *min = 0;
    *max = ln_prob.size() - 1;
  } else if (num_min == 1) {
    if (phase == 0) {
      *min = 0;
      *max = mins[0];
    } else if (phase == 1) {
      *min = mins[0];
      *max = ln_prob.size() - 1;
    } else {
      ERROR("unrecognized phase: " << phase);
    }
  } else {
    ERROR("multiple minima: " << num_min << " not implemented");
  }
}

bool FlatHistogram::is_equal(const FlatHistogram* flat_histogram,
    const double tolerance) const {
  if (!Criteria::is_equal(flat_histogram, tolerance)) {
      return false;
  }
  if (macrostate_old_ != flat_histogram->macrostate_old_) return false;
  if (macrostate_new_ != flat_histogram->macrostate_new_) return false;
  if (macrostate_current_ != flat_histogram->macrostate_current_) return false;
  return true;
}

}  // namespace feasst

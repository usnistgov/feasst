#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/histogram.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "system/include/system.h"
#include "monte_carlo/include/acceptance.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/wang_landau.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wltm.h"
#include "flat_histogram/include/macrostate_energy.h"
#include "flat_histogram/include/wang_landau.h"
#include "flat_histogram/include/bias.h"
#include "flat_histogram/include/macrostate.h"
#include "flat_histogram/include/flat_histogram.h"

namespace feasst {

FlatHistogram::FlatHistogram() {}  // do not use this constructor.
FlatHistogram::~FlatHistogram() {}

void FlatHistogram::init_(std::shared_ptr<Macrostate> macrostate,
    std::shared_ptr<Bias> bias) {
  macrostate_ = macrostate;
  bias_ = bias;
  bias_->resize(macrostate_->histogram());
}

FlatHistogram::FlatHistogram(std::shared_ptr<Macrostate> macrostate,
    std::shared_ptr<Bias> bias) {
  class_name_ = "FlatHistogram";
  init_(macrostate, bias);
}
FlatHistogram::FlatHistogram(argtype * args) : Criteria(args) {
  class_name_ = "FlatHistogram";
  ASSERT(!used("cycles_to_complete", *args),
    "FlatHistogram does not use the argument cycles_to_complete");
  init_(MacrostateEnergy().factory(str("Macrostate", args), args),
        MakeWangLandau({{"min_flatness", "1"}})->factory(str("Bias", args), args));
}
FlatHistogram::FlatHistogram(argtype args) : FlatHistogram(&args) {
  feasst_check_all_used(args);
}

FlatHistogram::FlatHistogram(std::shared_ptr<Macrostate> macrostate,
    std::shared_ptr<Bias> bias,
    std::shared_ptr<Constraint> constraint)
  : FlatHistogram(macrostate, bias) {
  add(constraint);
}

bool FlatHistogram::is_accepted(
    const System& system,
    Acceptance * acceptance,
    Random * random) {
  ASSERT(bias_ != NULL, "bias must be initialized before trials");
  bool is_accepted;
  double ln_metropolis_prob = acceptance->ln_metropolis_prob();
  DEBUG("macroshift " << acceptance->macrostate_shift());
//  const int shift = acceptance->macrostate_shift()*num_trial_states();
  if (acceptance->reject() ||
      !is_allowed(system, *acceptance) ) {
    is_accepted = false;
    ln_metropolis_prob = -NEAR_INFINITY;
    macrostate_new_ = macrostate_old_;
    DEBUG("forced rejection");
  } else {
    // the shift factory multiplied by number of states assumes only one
    // particle is added during an entire growth expanded cycle
    //macrostate_new_ = macrostate_->bin(system, this) + shift;
    macrostate_new_ = macrostate_->bin(system, *this, *acceptance);
    DEBUG("old " << macrostate_old_ << " new " << macrostate_new_);
    DEBUG("bias " << bias_->ln_bias(macrostate_new_, macrostate_old_));
    DEBUG("ln new " << bias_->ln_prob().value(macrostate_new_));
    DEBUG("ln old " << bias_->ln_prob().value(macrostate_old_));
    DEBUG("ln met " << ln_metropolis_prob);
    DEBUG("ln tot " << ln_metropolis_prob + bias_->ln_bias(macrostate_new_, macrostate_old_));
    if (macrostate_->is_allowed(system, *this, *acceptance) &&
        random->uniform() < exp(ln_metropolis_prob +
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
  bool is_endpoint = false;
  if (macrostate_new_ == macrostate_->soft_min() ||
      macrostate_new_ == macrostate_->soft_max()) {
    is_endpoint = true;
  }
  DEBUG(macrostate_old_ << " " <<  macrostate_new_ << " " << ln_metropolis_prob << " " << is_accepted << " " << is_endpoint);

  if (is_accepted) {
    ASSERT(system.num_configurations() == 1, "assumes 1 config");
    set_current_energy(acceptance->energy_new());
    set_current_energy_profile(acceptance->energy_profile_new());
    DEBUG("current energy: " << current_energy());
    macrostate_current_ = macrostate_new_;
  } else {
    // return the macrostate to the current value, as used by Analyze, etc.
    macrostate_current_ = macrostate_old_;
  }
  was_accepted_ = is_accepted;
  acceptance->set_endpoint(is_endpoint);
  return is_accepted;
}

std::string FlatHistogram::write() const {
  std::stringstream ss;
  ss << "#";
  ss << Criteria::write();
  ss << bias_->write();
  ss << "\"soft_min\":" << macrostate_->soft_min() << ","
     << "\"soft_max\":" << macrostate_->soft_max();
  ss << std::endl;
  ss << "state,"
     << bias_->write_per_bin_header("")
     << std::endl;
  const Histogram& hist = macrostate_->histogram();
  for (int bin = 0; bin < hist.size(); ++bin) {
//  for (int bin = macrostate_->soft_min();
//           bin <= macrostate_->soft_max();
//           ++bin) {
    ss << hist.center_of_bin(bin) << ","
       << bias_->write_per_bin(bin)
       << std::endl;
  }
  return ss.str();
}

void FlatHistogram::finalize(const Acceptance& acceptance) {
  DEBUG("macrostate_old_, " << macrostate_old_);
  DEBUG("macrostate_new_, " << macrostate_new_);
  DEBUG("acceptance.ln_metropolis_prob " << acceptance.ln_metropolis_prob());
  DEBUG("was_accepted_, " << was_accepted_);
  DEBUG("acceptance.endpoint " << acceptance.endpoint());
  bias_->update(macrostate_old_,
                macrostate_new_,
                acceptance.ln_metropolis_prob(),
                was_accepted_,
                acceptance.endpoint(),
                *macrostate_);
}

void FlatHistogram::revert_(const bool accepted, const bool endpoint, const double ln_prob, const std::vector<int>& updated) {
  Criteria::revert_(accepted, endpoint, ln_prob, updated);
//  if (!accepted) {
//    bias_->update_or_revert(macrostate_old_,
//                            macrostate_new_,
//                            ln_prob,
//                            accepted,
//                            endpoint,
//                            true);
//  }
  macrostate_new_ = macrostate_old_;
}

void FlatHistogram::imitate_trial_rejection_(const double ln_prob,
    const int state_old,
    const int state_new,
    const bool endpoint) {
  DEBUG("hi");
  DEBUG(state_old << " " <<  state_new << " " << ln_prob << " " << endpoint);
  bias_->update(state_old, state_new, ln_prob, false, endpoint, *macrostate_);
}

FEASST_MAPPER(FlatHistogram,);

FlatHistogram::FlatHistogram(std::istream& istr)
  : Criteria(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6079, "version mismatch: " << version);
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
  feasst_serialize_version(6079, ostr);
  feasst_serialize_fstdr(bias_, ostr);
  feasst_serialize_fstdr(macrostate_, ostr);
  feasst_serialize(macrostate_old_, ostr);
  feasst_serialize(macrostate_new_, ostr);
  feasst_serialize(macrostate_current_, ostr);
  feasst_serialize(is_macrostate_set_, ostr);
}

void FlatHistogram::before_attempt(const System& system) {
  if (!empty_) {
    empty_ = std::make_unique<Acceptance>();
  }
  macrostate_old_ = macrostate_->bin(system, *this, *empty_);
  DEBUG("macro old " << macrostate_old_);
  ASSERT(macrostate_old_ >= macrostate_->soft_min() &&
         macrostate_old_ <= macrostate_->soft_max(),
    "macrostate: " << macrostate_old_ << " is not in range from " <<
    macrostate_->soft_min() << " to " << macrostate_->soft_max());
}

bool FlatHistogram::is_fh_equal(const FlatHistogram& flat_histogram,
    const double tolerance) const {
  if (!Criteria::is_equal(flat_histogram, tolerance)) {
      return false;
  }
  if (macrostate_old_ != flat_histogram.macrostate_old_) return false;
  if (macrostate_new_ != flat_histogram.macrostate_new_) return false;
  if (macrostate_current_ != flat_histogram.macrostate_current_) return false;
  return true;
}

std::unique_ptr<FlatHistogram> FlatHistogram::flat_histogram(const Criteria& criteria) {
  std::stringstream ss;
  criteria.serialize(ss);
  return std::make_unique<FlatHistogram>(ss);
}

int FlatHistogram::set_soft_max(const int index, const System& sys) {
  if (!bias_->is_adjust_allowed(*macrostate_)) {
    return 0.;
  }
  return macrostate_->set_soft_max(index, sys, *this);
}

int FlatHistogram::set_soft_min(const int index, const System& sys) {
  if (!bias_->is_adjust_allowed(*macrostate_)) {
    return 0.;
  }
  return macrostate_->set_soft_min(index, sys, *this);
}

void FlatHistogram::set_cm(const bool inc_max, const int macro, const Criteria& crit) {
  if (inc_max) {
    macrostate_->add_to_soft_max(1);
  } else {
    macrostate_->remove_from_soft_min(1);
  }
  bias_->set_cm(macro, crit.bias());
}

void FlatHistogram::check_left_and_right_most_(const bool left_most, const bool right_most,
  const bool all_min_size,
  const int min_size, const System& system, const System * upper_sys,
  Criteria * criteria) {
  // HWH this only works if using new sweeps definition in TM
  if (left_most) {
    // consider increasing the soft_min if the existing visits are
    // enough for completion
    bool not_reject = true;
    while (not_reject) {
      not_reject = false;
      if (all_min_size && is_complete()) {
        set_soft_min(0, system);
      } else {
        if (bias().visits(macrostate().soft_min(), 0) >= cycles_to_complete() &&
            bias().visits(macrostate().soft_min(), 1) >= cycles_to_complete()) {
          if (macrostate().soft_max() - macrostate().soft_min() + 1 > min_size) {
            if (set_soft_min(macrostate().soft_min() + 1, system) != 0) {
              not_reject = true;
            }
          }
        }
      }
    }
  }
  if (right_most) {
    // consider decreasing the soft_max if the existing visits are
    // enough for completion
    bool not_reject = true;
    while (not_reject) {
      not_reject = false;
      if (all_min_size && criteria->is_complete()) {
        criteria->set_soft_max(criteria->num_states() - 1, *upper_sys);
      } else if (upper_sys) {
        if (criteria->bias().visits(criteria->macrostate().soft_max(), 0) >= criteria->cycles_to_complete() &&
            criteria->bias().visits(criteria->macrostate().soft_max(), 1) >= criteria->cycles_to_complete()) {
          if (criteria->macrostate().soft_max() - criteria->macrostate().soft_min() + 1 > min_size) {
            if (criteria->set_soft_max(criteria->macrostate().soft_max() - 1, *upper_sys) != 0) {
              not_reject = true;
            }
          }
        }
      } else {
        if (bias().visits(macrostate().soft_max(), 0) >= cycles_to_complete() &&
            bias().visits(macrostate().soft_max(), 1) >= cycles_to_complete()) {
          if (macrostate().soft_max() - macrostate().soft_min() + 1 > min_size) {
            if (set_soft_max(macrostate().soft_max() - 1, system) != 0) {
              not_reject = true;
            }
          }
        }
      }
    }
  }
}
void FlatHistogram::adjust_bounds(const bool left_most, const bool right_most,
  const bool left_complete, const bool right_complete,
  const bool all_min_size,
  const int min_size, const System& system, const System * upper_sys,
  Criteria * criteria, bool * adjusted_up, std::vector<int> * states) {
  DEBUG("left_most " << left_most);
  DEBUG("right_most " << right_most);
  check_left_and_right_most_(left_most, right_most, all_min_size, min_size, system, upper_sys, criteria);
  if (upper_sys) {
    *adjusted_up = false; // prevent hot potato adjustment. Only up, or only down, per adjust.
    bool not_reject = true;
    while (not_reject) {
      not_reject = false;
      const int lower_max = macrostate().soft_max();
      // if its left_most and already finished, don't send macrostates to upper
      if (num_cycles() < criteria->num_cycles()) {
//      && (!left_most || right_complete || all_min_size || num_cycles() < cycles_to_complete())) {
        if (lower_max - macrostate().soft_min() + 1 > min_size) {
          DEBUG("move macrostate from lower to upper");
          if (set_soft_max(lower_max - 1, system) > 0) {
            criteria->set_cm(false, lower_max, *this);
            not_reject = true;
            *adjusted_up = true;
            states->push_back(lower_max);
          }
        }
      }
      if (not_reject) {
        update();
        criteria->update();
      }
    }
    DEBUG("adjusted_up " << *adjusted_up);
    not_reject = true;
    while (not_reject && !*adjusted_up) {
      not_reject = false;
      const int upper_min = criteria->macrostate().soft_min();
      // if its right_most and already finished, don't send macrostates to lower
      if (num_cycles() > criteria->num_cycles()) {
//      && (!right_most || left_complete || all_min_size || criteria->num_cycles() < criteria->cycles_to_complete())) {
        if (criteria->macrostate().soft_max() - upper_min + 1 > min_size) {
          DEBUG("move macrostate from upper to lower");
          if (criteria->set_soft_min(upper_min + 1, *upper_sys) > 0) {
            set_cm(true, upper_min, *criteria);
            not_reject = true;
            states->push_back(upper_min);
          }
        }
      }
      if (not_reject) {
        update();
        criteria->update();
      }
    }
    DEBUG("states: " << feasst_str(*states));
  }
  check_left_and_right_most_(left_most, right_most, all_min_size, min_size, system, upper_sys, criteria);
}

int FlatHistogram::soft_min() const { return macrostate_->soft_min(); }
int FlatHistogram::soft_max() const { return macrostate_->soft_max(); }

void FlatHistogram::update_state(const System& system, const Acceptance& accept) {
  macrostate_current_ = macrostate_->bin(system, *this, accept);
}

int FlatHistogram::num_states() const { return macrostate_->histogram().size(); }
int FlatHistogram::phase() const { return bias_->phase(); }
void FlatHistogram::increment_phase() { bias_->increment_phase(); }
void FlatHistogram::set_ln_prob(const LnProbability& ln_prob) {
  bias_->set_ln_prob(ln_prob); }
const LnProbability& FlatHistogram::ln_prob() const { return bias_->ln_prob(); }
void FlatHistogram::update() { bias_->infrequent_update(*macrostate_); }
const Bias& FlatHistogram::bias() const { return const_cast<Bias&>(*bias_); }
void FlatHistogram::set_bias(std::shared_ptr<Bias> bias) { bias_ = bias; }
int FlatHistogram::cycles_to_complete() const {
  return bias_->cycles_to_complete(); }
void FlatHistogram::set_cycles_to_complete(const int num) {
  bias_->set_cycles_to_complete(num); }
int FlatHistogram::num_cycles(const int state) const {
  return bias_->num_cycles(state, *macrostate_); }
bool FlatHistogram::is_complete() const { return bias_->is_complete(); }
void FlatHistogram::set_complete() { bias_->set_complete_(); }
const Macrostate& FlatHistogram::macrostate() const {
  return const_cast<Macrostate&>(*macrostate_); }
int FlatHistogram::state() const { return macrostate_current_; }
int FlatHistogram::state_old() const { return macrostate_old_; }
int FlatHistogram::state_new() const { return macrostate_new_; }

}  // namespace feasst

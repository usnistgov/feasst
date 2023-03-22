
#include <algorithm>
#include "flat_histogram/include/wltm.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "utils/include/debug.h"

namespace feasst {

WLTM::WLTM(argtype * args) {
  class_name_ = "WLTM";
  collect_flatness_ = integer("collect_flatness", args);
  min_collect_sweeps_ = integer("min_collect_sweeps", args, -1);
  wang_landau_ = std::make_shared<WangLandau>(args);
  min_flatness_ = wang_landau_->min_flatness();
  ASSERT(collect_flatness_ < min_flatness_,
    "collect_flatness:" << collect_flatness_ << " should be less than " <<
    "transition_flatness:" << min_flatness_);
  transition_matrix_ = std::make_shared<TransitionMatrix>(args);
}
WLTM::WLTM(argtype args) : WLTM(&args) {
  FEASST_CHECK_ALL_USED(args);
}

bool WLTM::is_wl_bias_(const Macrostate& macro) const {
  if ((wang_landau_->num_flatness() < min_flatness_) ||
      (transition_matrix_->num_iterations(-1, macro) < min_collect_sweeps_)) {
    return true;
  }
  return false;
}

void WLTM::update(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool is_endpoint,
    const Macrostate& macro) {
  DEBUG("num_flatness " << wang_landau_->num_flatness());
  if (is_wl_bias_(macro)) {
    DEBUG("wl update");
    wang_landau_->update(macrostate_old, macrostate_new,
      ln_metropolis_prob, is_accepted, is_endpoint, macro);
  } else {
    if (production_ == 0) {
      production_ = 1;
      increment_phase();
    }
  }
  if (wang_landau_->num_flatness() >= collect_flatness_) {
    DEBUG("tm update");
    transition_matrix_->update(macrostate_old, macrostate_new,
      ln_metropolis_prob, is_accepted, is_endpoint, macro);
  }
  if (transition_matrix_->is_complete()) set_complete_();
}

const LnProbability& WLTM::ln_prob() const {
  if (wang_landau_->num_flatness() < min_flatness_) {
    return wang_landau_->ln_prob();
  } else {
    return transition_matrix_->ln_prob();
  }
}

void WLTM::resize(const Histogram& histogram) {
  wang_landau_->resize(histogram);
  transition_matrix_->resize(histogram);
}

void WLTM::infrequent_update(const Macrostate& macro) {
  if (is_wl_bias_(macro)) {
    wang_landau_->infrequent_update(macro);
  }
  if (wang_landau_->num_flatness() >= collect_flatness_) {
    transition_matrix_->infrequent_update(macro);
  }
}

std::string WLTM::write() const {
  std::stringstream ss;
  ss << Bias::write();
  ss << wang_landau_->write();
  ss << transition_matrix_->write();
  return ss.str();
}

std::string WLTM::write_per_bin_header() const {
  std::stringstream ss;
  ss << wang_landau_->write_per_bin_header() << ",";
  ss << transition_matrix_->write_per_bin_header();
  return ss.str();
}

std::string WLTM::write_per_bin(const int bin) const {
  std::stringstream ss;
  ss << wang_landau_->write_per_bin(bin) << ",";
  ss << transition_matrix_->write_per_bin(bin);
  return ss.str();
}

void WLTM::set_ln_prob(
    const LnProbability& ln_prob) {
  FATAL("not implemented. Try WangLandau or TransitionMatrix for reweighting.");
}

class MapWLTM {
 public:
  MapWLTM() {
    auto obj = MakeWLTM({{"collect_flatness", "1"},
                         {"min_flatness", "2"},
                         {"min_sweeps", "3"}});
    obj->deserialize_map()["WLTM"] = obj;
  }
};

static MapWLTM mapper_ = MapWLTM();

std::shared_ptr<Bias> WLTM::create(std::istream& istr) const {
  return std::make_shared<WLTM>(istr);
}

int WLTM::num_iterations(const int state, const Macrostate& macro) const {
  if (is_wl_bias_(macro)) {
    return 0;
  }
  return transition_matrix_->num_iterations(state, macro);
}

bool WLTM::is_adjust_allowed(const Macrostate& macro) const {
  if (is_wl_bias_(macro)) {
    return wang_landau_->is_adjust_allowed(macro);
  } else {
    return transition_matrix_->is_adjust_allowed(macro);
  }
}

WLTM::WLTM(std::istream& istr) : Bias(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1946, "mismatch version: " << version);
  feasst_deserialize(&collect_flatness_, istr);
  feasst_deserialize(&min_flatness_, istr);
  feasst_deserialize(&min_collect_sweeps_, istr);
  feasst_deserialize(&production_, istr);
  // HWH for unknown reasons, this function template does not work.
  // feasst_deserialize_fstdr(wang_landau_, istr);
  // feasst_deserialize_fstdr(transition_matrix_, istr);
  int existing;
  istr >> existing;
  if (existing != 0) {
    wang_landau_ = std::make_shared<WangLandau>(istr);
  }
  istr >> existing;
  if (existing != 0) {
    transition_matrix_ = std::make_shared<TransitionMatrix>(istr);
  }
}

void WLTM::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bias_(ostr);
  feasst_serialize_version(1946, ostr);
  feasst_serialize(collect_flatness_, ostr);
  feasst_serialize(min_flatness_, ostr);
  feasst_serialize(min_collect_sweeps_, ostr);
  feasst_serialize(production_, ostr);
  feasst_serialize_fstdr(wang_landau_, ostr);
  feasst_serialize_fstdr(transition_matrix_, ostr);
}

}  // namespace feasst


#include <algorithm>
#include "flat_histogram/include/wang_landau.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "utils/include/debug.h"
#include "flat_histogram/include/macrostate.h"

namespace feasst {

WangLandau::WangLandau(argtype args) : WangLandau(&args) {
  check_all_used(args);
}
WangLandau::WangLandau(argtype * args) {
  class_name_ = "WangLandau";
  min_flatness_ = integer("min_flatness", args);
  flatness_threshold_ = dble("flatness_threshold", args, 0.8);
  add_to_ln_probability_ = dble("add_to_ln_probability", args, 1.);
  reduce_ln_probability_ = dble("reduce_ln_probability", args, 0.5);
  min_visit_per_macro_ = integer("min_visit_per_macro", args, 1e3);
}

void WangLandau::flatness_update_() {
  DEBUG(feasst_str(visited_states_));
  std::fill(visited_states_.begin(), visited_states_.end(), 0);
  add_to_ln_probability_ *= reduce_ln_probability_;
  DEBUG("add_to_ln_probability_ " << add_to_ln_probability_);
  ++num_flatness_;
  ln_prob_.normalize();
  if (num_flatness_ >= min_flatness_) {
    set_complete_();
  }
}

void WangLandau::infrequent_update(const Macrostate& macro) {
  const int min_visit = *std::min_element(visited_states_.begin() + macro.soft_min(),
                                          visited_states_.begin() + macro.soft_max());
  const int num_visit = std::accumulate(visited_states_.begin() + macro.soft_min(),
                                        visited_states_.begin() + macro.soft_max(), 0);
  const double average = num_visit/static_cast<double>(macro.soft_max() - macro.soft_min() + 1);
  if ((min_visit >= min_visit_per_macro_) &&
      (min_visit >= flatness_threshold_ * average)) {
    flatness_update_();
  }
}

void WangLandau::update_or_revert(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool is_endpoint,
    const bool revert) {
  ASSERT(!revert, "not implemented");
  int bin = bin_(macrostate_old, macrostate_new, is_accepted);
  ln_prob_.add(bin, add_to_ln_probability_);
  ++visited_states_[bin];
}

void WangLandau::resize(const Histogram& histogram) {
  ln_prob_.resize(histogram.size());
  visited_states_.resize(histogram.size());
}

std::string WangLandau::write() const {
  std::stringstream ss;
  ss << Bias::write();
  ss << "num_flatness," << num_flatness_ << std::endl;
  return ss.str();
}

std::string WangLandau::write_per_bin_header() const {
  std::stringstream ss;
  ss << Bias::write_per_bin_header() << ",visited";
  return ss.str();
}

std::string WangLandau::write_per_bin(const int bin) const {
  std::stringstream ss;
  ss << Bias::write_per_bin(bin) << "," << visited_states_[bin];
  return ss.str();
}

void WangLandau::set_ln_prob(
    const LnProbability& ln_prob) {
  ASSERT(ln_prob.size() == ln_prob_.size(), "size mismatch: " <<
    ln_prob.size() << " " << ln_prob_.size());
  ln_prob_ = ln_prob;
}

class MapWangLandau {
 public:
  MapWangLandau() {
    auto obj = MakeWangLandau({{"min_flatness", "0"}});
    obj->deserialize_map()["WangLandau"] = obj;
  }
};

static MapWangLandau mapper_ = MapWangLandau();

WangLandau::WangLandau(std::istream& istr) : Bias(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 247, "mismatch version: " << version);
  feasst_deserialize_fstobj(&(ln_prob_), istr);
  feasst_deserialize(&(add_to_ln_probability_), istr);
  feasst_deserialize(&(reduce_ln_probability_), istr);
  feasst_deserialize(&(min_visit_per_macro_), istr);
  feasst_deserialize(&(flatness_threshold_), istr);
  feasst_deserialize(&(visited_states_), istr);
  feasst_deserialize(&(num_flatness_), istr);
  feasst_deserialize(&(min_flatness_), istr);
}

void WangLandau::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bias_(ostr);
  feasst_serialize_version(247, ostr);
  feasst_serialize_fstobj(ln_prob_, ostr);
  feasst_serialize(add_to_ln_probability_, ostr);
  feasst_serialize(reduce_ln_probability_, ostr);
  feasst_serialize(min_visit_per_macro_, ostr);
  feasst_serialize(flatness_threshold_, ostr);
  feasst_serialize(visited_states_, ostr);
  feasst_serialize(num_flatness_, ostr);
  feasst_serialize(min_flatness_, ostr);
}

void WangLandau::set_num_iterations_to_complete(const int flatness) {
  min_flatness_ = flatness;
  if (num_flatness_ < min_flatness_) set_incomplete_();
}

}  // namespace feasst

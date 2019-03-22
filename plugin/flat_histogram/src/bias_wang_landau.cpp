
#include <algorithm>
#include "flat_histogram/include/bias_wang_landau.h"
#include "core/include/utils_io.h"
#include "core/include/utils_math.h"
#include "core/include/debug.h"

namespace feasst {

BiasWangLandau::BiasWangLandau(const argtype &args) {
  args_.init(args);
  min_flatness_ = args_.key("min_flatness").integer();
  flatness_threshold_ = args_.key("flatness_threshold").dflt("0.8").dble();
  add_to_ln_probability_ = args_.key("add_to_ln_probability").dflt("1.").dble();
  reduce_ln_probability_ = args_.key("reduce_ln_probability").dflt("0.5").dble();
}

void BiasWangLandau::flatness_update_() {
  DEBUG(feasst_str(visited_states_));
  std::fill(visited_states_.begin(), visited_states_.end(), 0);
  add_to_ln_probability_ *= reduce_ln_probability_;
  ++num_flatness_;
  ln_macro_prob_.normalize();
  if (num_flatness_ >= min_flatness_) {
    set_complete_();
  }
}

void BiasWangLandau::update(const int macrostate_old,
                            const int macrostate_new,
                            const double ln_metropolis_prob,
                            const bool is_accepted) {
  int bin = bin_(macrostate_old, macrostate_new, is_accepted);
  ln_macro_prob_.add(bin, add_to_ln_probability_);
  ++visited_states_[bin];
  if (*std::min_element(visited_states_.begin(), visited_states_.end()) >
      flatness_threshold_ * average(visited_states_)) {
    flatness_update_();
  }
}

void BiasWangLandau::resize(const Histogram& histogram) {
  ln_macro_prob_.resize(histogram.size());
  visited_states_.resize(histogram.size());
}

std::string BiasWangLandau::write() const {
  std::stringstream ss;
  ss << Bias::write();
  ss << "num_flatness " << num_flatness_ << std::endl;
  return ss.str();
}

std::string BiasWangLandau::write_per_bin_header() const {
  std::stringstream ss;
  ss << Bias::write_per_bin_header() << " visited";
  return ss.str();
}

std::string BiasWangLandau::write_per_bin(const int bin) const {
  std::stringstream ss;
  ss << Bias::write_per_bin(bin) << " " << visited_states_[bin];
  return ss.str();
}

}  // namespace feasst

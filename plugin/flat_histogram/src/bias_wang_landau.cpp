
#include <algorithm>
#include "flat_histogram/include/bias_wang_landau.h"
#include "core/include/utils_io.h"

namespace feasst {

void BiasWangLandau::flatness_update_() {
  std::fill(visited_states_.begin(), visited_states_.end(), 0);
  ln_macro_prob_add_ *= ln_reduce_;
  ++flatness_;
}

void BiasWangLandau::update(const int macrostate_old_,
                            const int macrostate_new_,
                            const double ln_metropolis_prob_,
                            const bool is_accepted) {
  int bin = -1;
  if (is_accepted) {
    bin = macrostate_new_;
  } else {
    bin = macrostate_old_;
  }
  if (verbose == 1) cout << "bin " << bin << endl;
  ln_macro_prob_[bin] += ln_macro_prob_add_;
  ++visited_states_[bin];
  if (*std::min_element(visited_states_.begin(), visited_states_.end()) >
      visited_threshold_ * average(visited_states_)) {
    flatness_update_();
  }
}

void BiasWangLandau::resize(const Histogram& histogram) {
  ln_macro_prob_.resize(histogram.size());
  visited_states_.resize(histogram.size());
  if (verbose == 1) cout << "sizing bias " << histogram.size() << endl;
}

}  // namespace feasst

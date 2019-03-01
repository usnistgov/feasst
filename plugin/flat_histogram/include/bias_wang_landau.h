
#ifndef FEASST_FLAT_HISTOGRAM_BIAS_WANG_LANDAU_H_
#define FEASST_FLAT_HISTOGRAM_BIAS_WANG_LANDAU_H_

#include <vector>
#include <memory>
#include "flat_histogram/include/bias.h"
#include "core/include/histogram.h"
#include "core/include/utils_math.h"

namespace feasst {

/**
 Wang Landau flat histogram bias.
 HWH: Add literature reference
 */
class BiasWangLandau : public Bias {
 public:
  double ln_bias(const int bin_new, const int bin_old) const override {
    return ln_macro_prob_[bin_old] - ln_macro_prob_[bin_new];
  }

  void update(const int macrostate_old_,
              const int macrostate_new_,
              const double ln_metropolis_prob_,
              const bool is_accepted) override;

  std::vector<double> ln_macro_prob() const override { return ln_macro_prob_; }

  void resize(const Histogram& histogram);

  int verbose = 0;

  virtual ~BiasWangLandau() {}

 private:
  std::vector<double> ln_macro_prob_;

  /// Count of the number of times a state has been visited since the last time
  /// this histogram was reset after it was deemed to be sufficiently flat.
  std::vector<int> visited_states_;

  /// The amount to add to the natural log of the macrostate probability upon
  /// visiting that state.
  double ln_macro_prob_add_ = 1.;

  /// Reduce the amount to add to the natural log of the macrostate probability
  /// by this factor upon reaching a sufficiently flat histogram.
  const double ln_reduce_ = 0.5;

  /// The visited states histogram is determined to be flat when the percentage
  /// difference between minimum visisted states and average reaches this
  /// threshold.
  const double visited_threshold_ = 0.8;

  /// Number of times the visited states histogram was found to be flat.
  int flatness_ = 0;

  /// Perform update when the visited states histogram is found to be flat.
  void flatness_update_();
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_BIAS_WANG_LANDAU_H_

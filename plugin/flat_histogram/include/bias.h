
#ifndef FEASST_FLAT_HISTOGRAM_BIAS_H_
#define FEASST_FLAT_HISTOGRAM_BIAS_H_

namespace feasst {

/**
  Bias for flat histogram Monte Carlo.
  Assumes a one-dimensional macrostate.
 */
class Bias {
 public:
  /// Return the natural log of the bias for a transition from a macrostate
  /// in the old bin to a new bin.
  virtual double ln_bias(const int bin_new, const int bin_old) const = 0;

  /// Update the bias due to an attempted transition.
  virtual void update(const int macrostate_old_,
                      const int macrostate_new_,
                      const double ln_metropolis_prob_,
                      const bool is_accepted) = 0;

  /// The estimate of the natural log of the macrostate probability.
  virtual std::vector<double> ln_macro_prob() const = 0;

  virtual ~Bias() {}
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_BIAS_H_

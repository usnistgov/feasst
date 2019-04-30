
#ifndef FEASST_FLAT_HISTOGRAM_BIAS_H_
#define FEASST_FLAT_HISTOGRAM_BIAS_H_

#include <vector>
#include "core/include/histogram.h"
#include "core/include/random.h"
#include "flat_histogram/include/ln_probability_distribution.h"

namespace feasst {

/**
  Bias for flat histogram Monte Carlo.
  Assumes a one-dimensional macrostate.
 */
class Bias {
 public:
  Bias() {}

  /// Return the natural log of the bias for a transition from a macrostate
  /// in the old bin to a new bin.
  double ln_bias(const int bin_new, const int bin_old) const {
    return ln_macro_prob().value(bin_old) - ln_macro_prob().value(bin_new); }

  /// Update the bias due to an attempted transition.
  virtual void update(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted) = 0;

  /// The estimate of the natural log of the macrostate probability.
  virtual const LnProbabilityDistribution& ln_macro_prob() const = 0;

  virtual void resize(const Histogram& histogram) = 0;

  virtual std::string write() const { return std::string(""); }

  virtual std::string write_per_bin(const int bin) const;

  virtual std::string write_per_bin_header() const {
    return std::string("ln_prob"); }

  /// Return true if completion requirements are met.
  bool is_complete() const { return is_complete_; }

  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Bias> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Bias> >& deserialize_map();
  std::shared_ptr<Bias> deserialize(std::istream& istr);
  virtual ~Bias() {}

 protected:
  void set_complete_() { is_complete_ = true; }

  int bin_(
    const int macrostate_old,
    const int macrostate_new,
    const bool is_accepted);

  void serialize_bias_(std::ostream& ostr) const;
  Bias(std::istream& istr);

 private:
  bool is_complete_ = false;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_BIAS_H_

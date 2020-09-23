
#ifndef FEASST_FLAT_HISTOGRAM_BIAS_H_
#define FEASST_FLAT_HISTOGRAM_BIAS_H_

#include <vector>
#include "math/include/histogram.h"
#include "flat_histogram/include/ln_probability.h"

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
    return ln_prob().value(bin_old) - ln_prob().value(bin_new); }

  /// Update the bias due to an attempted transition.
  virtual void update_or_revert(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool revert) = 0;

  /// Update only.
  void update(
      const int macrostate_old,
      const int macrostate_new,
      const double ln_metropolis_prob,
      const bool is_accepted) {
    update_or_revert(macrostate_old, macrostate_new,
           ln_metropolis_prob, is_accepted, false);
  }


  /// Perform an infrequent update to the bias.
  virtual void infrequent_update() {}

  /// The natural log of the macrostate probability.
  virtual const LnProbability& ln_prob() const = 0;

  virtual void resize(const Histogram& histogram) = 0;

  /// Set the macrostate probability distribution.
  virtual void set_ln_prob(const LnProbability& ln_prob) = 0;

  virtual std::string write() const { return std::string(""); }

  virtual std::string write_per_bin(const int bin) const;

  virtual std::string write_per_bin_header() const {
    return std::string("ln_prob"); }

  /// Return true if completion requirements are met.
  bool is_complete() const { return is_complete_; }

  /// Return the simulation phase index used to differentiate production
  /// and initialization, etc.
  virtual int phase() const { return phase_; }

  /// Increment the simulation phase.
  virtual void increment_phase() { ++phase_; }

  /// Set the number of iterations required for completion.
  /// In TransitionMatrix and WLTM, this is the minimum number of sweeps.
  /// In WangLandau, this is the minimum number of flatness checks.
  /// Afterward, check again for completeness.
  virtual void set_num_iterations(const int iteration) = 0;

  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Bias> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Bias> >& deserialize_map();
  std::shared_ptr<Bias> deserialize(std::istream& istr);
  virtual ~Bias() {}

 protected:
  std::string class_name_ = "Bias";

  void set_complete_() { is_complete_ = true; }
  void set_incomplete_() { is_complete_ = false; }

  int bin_(
    const int macrostate_old,
    const int macrostate_new,
    const bool is_accepted);

  void serialize_bias_(std::ostream& ostr) const;
  Bias(std::istream& istr);

 private:
  bool is_complete_ = false;
  int phase_ = 0;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_BIAS_H_

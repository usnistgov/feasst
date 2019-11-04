
#ifndef FEASST_FLAT_HISTOGRAM_FLAT_HISTOGRAM_H_
#define FEASST_FLAT_HISTOGRAM_FLAT_HISTOGRAM_H_

#include <memory>
#include "math/include/accumulator.h"
#include "monte_carlo/include/criteria.h"
#include "flat_histogram/include/macrostate.h"
#include "flat_histogram/include/bias.h"

namespace feasst {

/**
  The macrostate must be defined before the bias.
  Use MacrostateAccumulator to compute custom per-macrostate quantities.
 */
class FlatHistogram : public Criteria {
 public:
  FlatHistogram(const argtype &args = argtype()) : Criteria(args) {}

  void before_attempt(const System* system) override {
    macrostate_old_ = macrostate_->bin(system, this);
    DEBUG("macro old " << macrostate_old_);
  }

  bool is_accepted(const Acceptance& acceptance,
    const System * system,
    const double uniform_random) override;

  std::string write() const override;
  bool is_complete() override { return bias_->is_complete(); }

  /// Set the macrostate which is subject to the bias.
  void set(const std::shared_ptr<Macrostate> macrostate) {
    macrostate_ = macrostate;
    is_macrostate_set_ = true;
  }

  /// Return the macrostate.
  const Macrostate * macrostate() const { return macrostate_.get(); }

  /// Set the bias for the flat histogram method.
  void set(const std::shared_ptr<Bias> bias) {
    ASSERT(is_macrostate_set_, "set macrostate before bias");
    bias_ = bias;
    bias_->resize(macrostate_->histogram());
  }

  /// Return the bias.
  const Bias * bias() const { return bias_.get(); }

  /// Return the state. Return -1 if state is not determined.
  int state() const override { return macrostate_new_; }
  int num_states() const override { return macrostate_->histogram().size(); }

  /// Revert changes from previous trial.
  void revert(const bool accepted) override;

  // HWH consider moving the below functions to a new class or util

  /// Return the pressure.
  /// Return a reweighted macrostate distribution.
  LnProbability reweight(
      /**
        Change in conjugate variable. For example, if reweighting in number
        of particles using the grand canonical ensemble, the change in the
        thermodynamic conjugate is \f$\Delta(\beta\mu\)f$.
        If reweighting in potential energy in the microcanonical, the change
        in the thermodynamic conjugate is \f$\Delta(-\beta)\f$.
       */
      const double delta_conjugate) {
    LnProbability lnpirw = deep_copy(bias()->ln_prob());
    for (int macro = 0; macro < lnpirw.size(); ++macro) {
      lnpirw.add(macro, macrostate()->histogram().center_of_bin(macro)
                 *delta_conjugate);
    }
    lnpirw.normalize();
    return lnpirw;
  }

  /// Set the macrostate probability distribution.
  void set_ln_prob(const LnProbability& ln_prob) {
    bias_->set_ln_prob(ln_prob); }

  // HWH add reference to Shen/Errington
  double pressure(const double volume,
      /// Select phase by order of macrostate.
      /// Assumes default method of dividing phase boundary.
      const int phase = 0) const {
    int min, max;
    phase_boundary_(phase, &min, &max);
    return (-ln_prob_().value(0) + log(ln_prob_().sum_probability(min, max)))
           /volume/beta();
  }

  /// Return the ensemble averaged property from a list of properties averaged
  /// at each macrostate.
  double average(const LnProbability& ln_prob,
                 // HWH make macrostate_averages work with bin averages.
                 const std::vector<double>& macrostate_averages,
                 /// Select phase by order of macrostate.
                 /// Assumes default method of dividing phase boundary.
                 const int phase = 0) const {
    ASSERT(ln_prob.size() == static_cast<int>(macrostate_averages.size()),
      "size mismatch: ln_prob:" << ln_prob.size() <<
      " macro:" << macrostate_averages.size());
    int min, max;
    phase_boundary_(ln_prob, phase, &min, &max);
    double average = 0.;
    for (int bin = min; bin < max + 1; ++bin) {
      average += macrostate_averages[bin]*exp(ln_prob.value(bin));
    }
    return average/ln_prob.sum_probability(min, max);
  }

  double average(const std::vector<double>& macrostate_averages,
                 /// Select phase by order of macrostate.
                 /// Assumes default method of dividing phase boundary.
                 const int phase = 0) const {
  return average(ln_prob_(), macrostate_averages, phase); }

  /// Return the average macrostate
  double average_macrostate(const LnProbability& ln_prob,
      /// Select phase by order of macrostate.
      /// Assumes default method of dividing phase boundary.
      const int phase = 0) const {
    int min, max;
    phase_boundary_(phase, &min, &max);
    double average = 0.;
    for (int bin = min; bin < max + 1; ++bin) {
      average += macrostate_->value(bin)*exp(ln_prob.value(bin));
    }
    return average/ln_prob.sum_probability(min, max);
  }

  double average_macrostate(
      /// Select phase by order of macrostate.
      /// Assumes default method of dividing phase boundary.
      const int phase = 0) const {
    return average_macrostate(ln_prob_(), phase);
  }

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<FlatHistogram>(istr); }
  void serialize(std::ostream& ostr) const override;
  FlatHistogram(std::istream& istr);
  ~FlatHistogram() {}

 private:
  const std::string class_name_ = "FlatHistogram";
  std::shared_ptr<Bias> bias_;
  std::shared_ptr<Macrostate> macrostate_;
  int macrostate_old_ = -1;
  int macrostate_new_ = -1;
  bool is_macrostate_set_ = false;

  const LnProbability& ln_prob_() const { return bias_->ln_prob(); }

  /// Determine min and max indices for a given phase
  /// Return -1 if no phase boundary.
  void phase_boundary_(const LnProbability& ln_prob,
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

  /// Same as above but with the ln_prob_ contained in this class.
  void phase_boundary_(const int phase, int * min, int * max) const {
    phase_boundary_(ln_prob_(), phase, min, max);
  }
};

inline std::shared_ptr<FlatHistogram> MakeFlatHistogram(
    const argtype &args = argtype()) {
  return std::make_shared<FlatHistogram>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_FLAT_HISTOGRAM_H_

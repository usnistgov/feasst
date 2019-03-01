
#ifndef FEASST_FLAT_HISTOGRAM_CRITERIA_FLAT_HISTOGRAM_H_
#define FEASST_FLAT_HISTOGRAM_CRITERIA_FLAT_HISTOGRAM_H_

#include <memory>
#include "core/include/criteria.h"
#include "flat_histogram/include/macrostate.h"
#include "flat_histogram/include/bias.h"

namespace feasst {

/**
 */
class CriteriaFlatHistogram : public Criteria {
 public:
  void before_attempt(const System* system) override {
    macrostate_old_ = macrostate_->bin(system, this);
    if (verbose) std::cout << "macro old " << macrostate_old_ << std::endl;
  }

  bool is_accepted(const AcceptanceCriteria accept_criteria) override;

  /// Set the bias for the flat histogram method.
  void set_bias(const std::shared_ptr<Bias> bias) { bias_ = bias; }

  /// Set the macrostate which is subject to the bias.
  void set_macrostate(const std::shared_ptr<Macrostate> macrostate) {
    macrostate_ = macrostate;
  }

  bool verbose = false;

 private:
  std::shared_ptr<Bias> bias_;
  std::shared_ptr<Macrostate> macrostate_;
  int macrostate_old_, macrostate_new_;
  Random random_;

  void after_attempt_(const System* system) override {
    macrostate_new_ = macrostate_->bin(system, this);
  }
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_CRITERIA_FLAT_HISTOGRAM_H_

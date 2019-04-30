
#ifndef FEASST_FLAT_HISTOGRAM_CRITERIA_FLAT_HISTOGRAM_H_
#define FEASST_FLAT_HISTOGRAM_CRITERIA_FLAT_HISTOGRAM_H_

#include <memory>
#include "core/include/accumulator.h"
#include "core/include/criteria.h"
#include "flat_histogram/include/macrostate.h"
#include "flat_histogram/include/bias.h"
#include "flat_histogram/include/macrostate_accumulator.h"

namespace feasst {

/**
  The macrostate must be defined before the bias.
  Use MacrostateAccumulator to compute custom per-macrostate quantities.
 */
class CriteriaFlatHistogram : public Criteria {
 public:
  CriteriaFlatHistogram(const argtype &args = argtype()) : Criteria(args) {
    bin_trackers_.add(std::make_shared<BinEnergy>());
  }

  void before_attempt(const System* system) override {
    macrostate_old_ = macrostate_->bin(system, this);
    DEBUG("macro old " << macrostate_old_);
    bin_trackers_.update(macrostate_old_, system, this);
  }

  bool is_accepted(const AcceptanceCriteria accept_criteria) override;
  std::string write() const override;
  bool is_complete() override { return bias_->is_complete(); }

  /// Set the bias for the flat histogram method.
  void set(const std::shared_ptr<Bias> bias) {
    ASSERT(is_macrostate_set_, "set macrostate before bias");
    bias_ = bias;
    bias_->resize(macrostate_->histogram());
    bin_trackers_.resize(macrostate_->histogram().size());
  }

  /// Set the macrostate which is subject to the bias.
  void set(const std::shared_ptr<Macrostate> macrostate) {
    macrostate_ = macrostate;
    is_macrostate_set_ = true;
  }

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<CriteriaFlatHistogram>(istr); }

  void serialize(std::ostream& ostr) const override;
  CriteriaFlatHistogram(std::istream& istr);
  ~CriteriaFlatHistogram() {}

 private:
  const std::string class_name_ = "CriteriaFlatHistogram";
  std::shared_ptr<Bias> bias_;
  std::shared_ptr<Macrostate> macrostate_;
  int macrostate_old_, macrostate_new_;
  MacrostateAccumulatorFactory bin_trackers_;
  bool is_macrostate_set_ = false;

  Random random_;

  void after_attempt_(const System* system) override {
    macrostate_new_ = macrostate_->bin(system, this);
  }
};

inline std::shared_ptr<CriteriaFlatHistogram> MakeCriteriaFlatHistogram(
    const argtype &args = argtype()) {
  return std::make_shared<CriteriaFlatHistogram>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_CRITERIA_FLAT_HISTOGRAM_H_


#ifndef FEASST_CHAIN_ANALYZE_BONDS_H_
#define FEASST_CHAIN_ANALYZE_BONDS_H_

#include <vector>
#include "math/include/accumulator.h"
#include "math/include/histogram.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Compute the distribution and moments of each type of bond and angle.
 */
class AnalyzeBonds : public AnalyzeUpdateOnly {
 public:
  /**
    args:
    - angle_bin_width: bin width for angle histogram in units of degrees.
      (default: 1).
    - bond_bin_width: bin width for bond histogram (default: 1).
   */
  AnalyzeBonds(const argtype &args = argtype());

  void update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) override;
  void initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;

  /// Return the average bond in degrees.
  const Accumulator& bond(const int type) const { return bond_[type]; }

  /// Return the average angle in degrees.
  const Accumulator& angle(const int type) const { return angle_[type]; }

  /// Return the histogram of bond in degrees.
  const Histogram& bond_hist(const int type) const { return bond_hist_[type]; }

  /// Return the histogram of angle in degrees.
  const Histogram& angle_hist(const int type) const {
    return angle_hist_[type]; }

  std::string class_name() const override { return std::string("AnalyzeBonds"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<AnalyzeBonds>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit AnalyzeBonds(std::istream& istr);

 private:
  std::vector<Accumulator> bond_, angle_;
  std::vector<Histogram> bond_hist_, angle_hist_;
};

inline std::shared_ptr<AnalyzeBonds> MakeAnalyzeBonds(const argtype &args = argtype()) {
  return std::make_shared<AnalyzeBonds>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_ANALYZE_BONDS_H_

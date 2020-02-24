
#ifndef FEASST_CHAIN_COMPUTE_DEPROTONATION_H_
#define FEASST_CHAIN_COMPUTE_DEPROTONATION_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

class ComputeDeprotonation : public TrialCompute {
 public:
  /**
    args:
    - pKa: the negative log-base-10 for the deprotonation reaction:
      \f$K_a = \frac{a_{XO}a_H}{a_{XOH}} \f$
      The default value, 0, assumes the input pH in Criteria already takes the
      pKa into account (e.g., pH was entered as the value of pH - pKa).
   */
  ComputeDeprotonation(const argtype& args = argtype());

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeDeprotonation(std::istream& istr);
  virtual ~ComputeDeprotonation() {}

 protected:
  void serialize_compute_deprotonation_(std::ostream& ostr) const;

 private:
  std::vector<TrialStage*> stage0;
  double pKa_;
};

inline std::shared_ptr<ComputeDeprotonation> MakeComputeDeprotonation(
    const argtype& args = argtype()) {
  return std::make_shared<ComputeDeprotonation>(args);
}
}  // namespace feasst

#endif  // FEASST_CHAIN_COMPUTE_DEPROTONATION_H_

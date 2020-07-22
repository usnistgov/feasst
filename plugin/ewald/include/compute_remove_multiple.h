
#ifndef FEASST_MONTE_CARLO_COMPUTE_REMOVE_MULTIPLE_H_
#define FEASST_MONTE_CARLO_COMPUTE_REMOVE_MULTIPLE_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

// HWH do away with macrostate shifts.. make easy add/delete.. finalize optim
// HWH parse args between classes... warn if not used.
/**
  Attempt to remove multiple particles.
  For a derivation of the acceptance criteria, see ComputeAddMultiple that is
  the reverse of this trial.
 */
class ComputeRemoveMultiple : public TrialCompute {
 public:
  /**
    args:
    - shift: macrostate shift (default: -1).
   */
  ComputeRemoveMultiple(const argtype& args = argtype());

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeRemoveMultiple(std::istream& istr);
  virtual ~ComputeRemoveMultiple() {}

 protected:
  void serialize_compute_remove_multiple_(std::ostream& ostr) const;
  int shift_;

  // temporary, not serialized
  std::vector<int> delta_;
};

inline std::shared_ptr<ComputeRemoveMultiple> MakeComputeRemoveMultiple() {
  return std::make_shared<ComputeRemoveMultiple>();
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_COMPUTE_REMOVE_MULTIPLE_H_

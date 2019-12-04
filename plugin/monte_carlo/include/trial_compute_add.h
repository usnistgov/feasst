
#ifndef FEASST_MONTE_CARLO_TRIAL_COMPUTE_ADD_H_
#define FEASST_MONTE_CARLO_TRIAL_COMPUTE_ADD_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

class TrialComputeAdd : public TrialCompute {
 public:
  TrialComputeAdd();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeAdd(std::istream& istr);
  virtual ~TrialComputeAdd() {}

 protected:
  void serialize_trial_compute_add_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialComputeAdd> MakeTrialComputeAdd() {
  return std::make_shared<TrialComputeAdd>();
}
}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_COMPUTE_ADD_H_

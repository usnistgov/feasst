
#ifndef FEASST_MONTE_CARLO_TRIAL_COMPUTE_ADD_REMOVE_H_
#define FEASST_MONTE_CARLO_TRIAL_COMPUTE_ADD_REMOVE_H_

#include <memory>
#include <vector>
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

class TrialComputeAdd;
class TrialComputeRemove;

/**
 */
class TrialComputeAddRemove : public TrialCompute {
 public:
  explicit TrialComputeAddRemove(argtype args = argtype());
  explicit TrialComputeAddRemove(argtype * args);
  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeAddRemove(std::istream& istr);
  virtual ~TrialComputeAddRemove();

 private:
  std::unique_ptr<TrialComputeAdd> add_;
  std::unique_ptr<TrialComputeRemove> rm_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_COMPUTE_ADD_REMOVE_H_

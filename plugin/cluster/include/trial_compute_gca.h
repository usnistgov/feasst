
#ifndef FEASST_CLUSTER_TRIAL_COMPUTE_GCA_H_
#define FEASST_CLUSTER_TRIAL_COMPUTE_GCA_H_

#include <vector>
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {
/*
  GCA needs a few novel bits
  - Make this Trial, not TrialCompute
  - Acceptance/Criteria need a force acceptance option?
  - how to recursive pivot while maintaining energymap?
    - energy change is only calculated from refused pivots
    - any accepted pivot would be finalized immediately before
      recursive pivot again
    - but before finalize, those particles affected by pivot would be
      added to the pivot list
 */
class TrialComputeGCA : public TrialCompute {
 public:
  TrialComputeGCA();

  void perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) override;
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeGCA(std::istream& istr);
  virtual ~TrialComputeGCA() {}
};

inline std::shared_ptr<TrialComputeGCA> MakeTrialComputeGCA() {
  return std::make_shared<TrialComputeGCA>();
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_COMPUTE_GCA_H_

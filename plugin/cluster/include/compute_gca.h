
#ifndef FEASST_CLUSTER_COMPUTE_GCA_H_
#define FEASST_CLUSTER_COMPUTE_GCA_H_

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
class ComputeGCA : public TrialCompute {
 public:
  ComputeGCA();

  void perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) override;
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeGCA(std::istream& istr);
  virtual ~ComputeGCA() {}
};

inline std::shared_ptr<ComputeGCA> MakeComputeGCA() {
  return std::make_shared<ComputeGCA>();
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_COMPUTE_GCA_H_


#ifndef FEASST_MONTE_CARLO_TRIAL_COMPUTE_REMOVE_H_
#define FEASST_MONTE_CARLO_TRIAL_COMPUTE_REMOVE_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

/**
  Attempt to remove a particle.
  See TrialComputeAdd for derivation of the acceptance probability that is
  the reverse of this Trial.
 */
class TrialComputeRemove : public TrialCompute {
 public:
  TrialComputeRemove();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeRemove(std::istream& istr);
  virtual ~TrialComputeRemove() {}

 protected:
  void serialize_trial_compute_remove_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialComputeRemove> MakeTrialComputeRemove() {
  return std::make_shared<TrialComputeRemove>();
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_COMPUTE_REMOVE_H_


#ifndef FEASST_MONTE_CARLO_TRIAL_COMPUTE_MOVE_CLUSTER_H_
#define FEASST_MONTE_CARLO_TRIAL_COMPUTE_MOVE_CLUSTER_H_

#include <vector>
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"
#include "cluster/include/trial_select_cluster.h"

namespace feasst {

class TrialComputeMoveCluster : public TrialCompute {
 public:
  TrialComputeMoveCluster();

  void perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) override;
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeMoveCluster(std::istream& istr);
  virtual ~TrialComputeMoveCluster() {}

 protected:
  void serialize_trial_compute_move_cluster_(std::ostream& ostr) const;

  // temporary
  std::shared_ptr<TrialSelectCluster> cselect_;
};

inline std::shared_ptr<TrialComputeMoveCluster> MakeTrialComputeMoveCluster() {
  return std::make_shared<TrialComputeMoveCluster>();
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_COMPUTE_MOVE_CLUSTER_H_

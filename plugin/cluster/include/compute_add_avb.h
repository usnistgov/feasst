
#ifndef FEASST_CLUSTER_COMPUTE_ADD_AVB_H_
#define FEASST_CLUSTER_COMPUTE_ADD_AVB_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

class ComputeAddAVB : public TrialCompute {
 public:
  ComputeAddAVB();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeAddAVB(std::istream& istr);
  virtual ~ComputeAddAVB() {}

 protected:
  void serialize_compute_add_avb_(std::ostream& ostr) const;
};

inline std::shared_ptr<ComputeAddAVB> MakeComputeAddAVB() {
  return std::make_shared<ComputeAddAVB>();
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_COMPUTE_ADD_AVB_H_

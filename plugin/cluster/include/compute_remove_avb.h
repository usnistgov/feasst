
#ifndef FEASST_CLUSTER_COMPUTE_REMOVE_AVB_H_
#define FEASST_CLUSTER_COMPUTE_REMOVE_AVB_H_

#include <memory>
#include <vector>
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class ComputeRemoveAVB : public TrialCompute {
 public:
  ComputeRemoveAVB();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeRemoveAVB(std::istream& istr);
  virtual ~ComputeRemoveAVB() {}

 protected:
  void serialize_compute_remove_avb_(std::ostream& ostr) const;
};

inline std::shared_ptr<ComputeRemoveAVB> MakeComputeRemoveAVB() {
  return std::make_shared<ComputeRemoveAVB>();
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_COMPUTE_REMOVE_AVB_H_

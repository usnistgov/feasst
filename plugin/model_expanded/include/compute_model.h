
#ifndef FEASST_MODEL_EXPANDED_COMPUTE_MODEL_H_
#define FEASST_MODEL_EXPANDED_COMPUTE_MODEL_H_

#include <memory>
#include <vector>
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class ComputeModel : public TrialCompute {
 public:
  ComputeModel();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeModel(std::istream& istr);
  virtual ~ComputeModel() {}

 protected:
  void serialize_trial_compute_add_(std::ostream& ostr) const;
};

inline std::shared_ptr<ComputeModel> MakeComputeModel() {
  return std::make_shared<ComputeModel>();
}
}  // namespace feasst

#endif  // FEASST_MODEL_EXPANDED_COMPUTE_MODEL_H_

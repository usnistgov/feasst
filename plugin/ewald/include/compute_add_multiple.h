
#ifndef FEASST_MONTE_CARLO_COMPUTE_ADD_MULTIPLE_H_
#define FEASST_MONTE_CARLO_COMPUTE_ADD_MULTIPLE_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

class ComputeAddMultiple : public TrialCompute {
 public:
  ComputeAddMultiple();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeAddMultiple(std::istream& istr);
  virtual ~ComputeAddMultiple() {}

 protected:
  void serialize_compute_add_multiple_(std::ostream& ostr) const;
};

inline std::shared_ptr<ComputeAddMultiple> MakeComputeAddMultiple() {
  return std::make_shared<ComputeAddMultiple>();
}
}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_COMPUTE_ADD_MULTIPLE_H_

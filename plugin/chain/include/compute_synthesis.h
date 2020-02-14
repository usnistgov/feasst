
#ifndef FEASST_CHAIN_COMPUTE_SYNTHESIS_H_
#define FEASST_CHAIN_COMPUTE_SYNTHESIS_H_

#include <vector>
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

class ComputeSynthesis : public TrialCompute {
 public:
  ComputeSynthesis();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeSynthesis(std::istream& istr);
  virtual ~ComputeSynthesis() {}

 protected:
  void serialize_compute_synthesis_(std::ostream& ostr) const;
};

inline std::shared_ptr<ComputeSynthesis> MakeComputeSynthesis() {
  return std::make_shared<ComputeSynthesis>();
}
}  // namespace feasst

#endif  // FEASST_CHAIN_COMPUTE_SYNTHESIS_H_

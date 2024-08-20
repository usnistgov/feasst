
#ifndef FEASST_MONTE_CARLO_COMPUTE_AVB4_H_
#define FEASST_MONTE_CARLO_COMPUTE_AVB4_H_

#include <memory>
#include <vector>
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute_move.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
 */
class ComputeAVB4 : public TrialComputeMove {
 public:
  ComputeAVB4();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeAVB4(std::istream& istr);
  virtual ~ComputeAVB4() {}

 protected:
  void serialize_compute_avb4_(std::ostream& ostr) const;
};

inline std::shared_ptr<ComputeAVB4> MakeComputeAVB4() {
  return std::make_shared<ComputeAVB4>();
}
}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_COMPUTE_AVB4_H_

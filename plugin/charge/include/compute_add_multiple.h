
#ifndef FEASST_CHARGE_COMPUTE_ADD_MULTIPLE_H_
#define FEASST_CHARGE_COMPUTE_ADD_MULTIPLE_H_

#include <memory>
#include <vector>
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class ComputeAddMultiple : public TrialCompute {
 public:
  /**
    args:
    - shift: macrostate shift (default: -1).
   */
  explicit ComputeAddMultiple(argtype args = argtype());
  explicit ComputeAddMultiple(argtype * args);

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
  int shift_;

 private:
  std::vector<int> delta_;
};

inline std::shared_ptr<ComputeAddMultiple> MakeComputeAddMultiple() {
  return std::make_shared<ComputeAddMultiple>();
}
}  // namespace feasst

#endif  // FEASST_CHARGE_COMPUTE_ADD_MULTIPLE_H_

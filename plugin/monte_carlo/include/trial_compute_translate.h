
#ifndef FEASST_MONTE_CARLO_TRIAL_COMPUTE_TRANSLATE_H_
#define FEASST_MONTE_CARLO_TRIAL_COMPUTE_TRANSLATE_H_

#include <vector>
#include "system/include/system.h"
#include "configuration/include/select.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

/**
  Translate a selection of particles and sites.
 */
class TrialComputeTranslate : public TrialCompute {
 public:
  explicit TrialComputeTranslate(argtype args = argtype());
  explicit TrialComputeTranslate(argtype * args);

  void perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) override;
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeTranslate(std::istream& istr);
  virtual ~TrialComputeTranslate() {}

 protected:
  void serialize_trial_compute_translate_(std::ostream& ostr) const;

 private:
  Select new_;
};

inline std::shared_ptr<TrialComputeTranslate> MakeTrialComputeTranslate(
    argtype args = argtype()) {
  return std::make_shared<TrialComputeTranslate>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_COMPUTE_TRANSLATE_H_

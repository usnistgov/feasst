
#ifndef FEASST_MONTE_CARLO_TRIAL_COMPUTE_H_
#define FEASST_MONTE_CARLO_TRIAL_COMPUTE_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

/// Implement the perturbation and calculation of acceptance.
class TrialCompute {
 public:
  TrialCompute() {}

  /// Perform the stages on the system and compute the acceptance.
  void compute_rosenbluth(
    /// Set to 1 for "old" system and "0" for new.
    const int old,
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random);

  /// Perform the Perturbations and determine acceptance.
  virtual void perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) = 0;

  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<TrialCompute> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<TrialCompute> >& deserialize_map();
  std::shared_ptr<TrialCompute> deserialize(std::istream& istr);
  virtual ~TrialCompute() {}

 protected:
  std::string class_name_ = "TrialCompute";

  void serialize_trial_compute_(std::ostream& ostr) const;
  TrialCompute(std::istream& istr);
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_COMPUTE_H_

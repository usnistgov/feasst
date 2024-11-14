
#ifndef FEASST_BETA_EXPANDED_COMPUTE_BETA_H_
#define FEASST_BETA_EXPANDED_COMPUTE_BETA_H_

#include <memory>
#include <vector>
#include <map>
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

class Criteria;
class System;
class TrialStage;

typedef std::map<std::string, std::string> argtype;

class ComputeBeta : public TrialCompute {
 public:
  ComputeBeta();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeBeta(std::istream& istr);
  virtual ~ComputeBeta() {}

 protected:
  void serialize_trial_compute_add_(std::ostream& ostr) const;
};

}  // namespace feasst

#endif  // FEASST_BETA_EXPANDED_COMPUTE_BETA_H_

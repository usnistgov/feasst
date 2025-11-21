
#ifndef FEASST_GIBBS_COMPUTE_GIBBS_MORPH_H_
#define FEASST_GIBBS_COMPUTE_GIBBS_MORPH_H_

#include <memory>
#include <vector>
#include <map>
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

class Criteria;
class System;
class TrialStage;

typedef std::map<std::string, std::string> argtype;

/**
  See TrialGibbsMorph.
 */
class ComputeGibbsMorph : public TrialCompute {
 public:
  ComputeGibbsMorph();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeGibbsMorph(std::istream& istr);
  virtual ~ComputeGibbsMorph() {}

 protected:
  void serialize_compute_gibbs_morph_(std::ostream& ostr) const;
};

}  // namespace feasst

#endif  // FEASST_GIBBS_COMPUTE_GIBBS_MORPH_H_

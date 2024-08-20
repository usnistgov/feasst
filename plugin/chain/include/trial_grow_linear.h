
#ifndef FEASST_CHAIN_TRIAL_GROW_LINEAR_H_
#define FEASST_CHAIN_TRIAL_GROW_LINEAR_H_

#include <memory>
#include "monte_carlo/include/trial.h"

namespace feasst {

/**
  Grow a freely jointed linear chain in multiple stages.
 */
class TrialGrowLinear : public Trial {
 public:
  /**
    args:
    - particle_type: type of particle in configuration (default: 0).
   */
  TrialGrowLinear(
    std::shared_ptr<TrialCompute> compute,
    argtype args = argtype());
  void precompute(Criteria * criteria, System * system) override;
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialGrowLinear(std::istream& istr);
  virtual ~TrialGrowLinear() {}

 private:
  argtype stored_args_;
};

inline std::shared_ptr<TrialGrowLinear> MakeTrialGrowLinear(
    std::shared_ptr<TrialCompute> compute,
    argtype args = argtype()) {
  return std::make_shared<TrialGrowLinear>(compute, args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_GROW_LINEAR_H_

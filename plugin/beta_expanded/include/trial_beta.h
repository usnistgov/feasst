#ifndef FEASST_BETA_EXPANDED_TRIAL_BETA_H_
#define FEASST_BETA_EXPANDED_TRIAL_BETA_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/**
  Attempt to change the inverse temperature, \f$\beta=\frac{1}{k_B T}\f$ by a
  fixed amount.
 */
class TrialBeta : public Trial {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - Trial arguments.
    - PerturbBeta arguments.
    - Tunable arguments.
   */
  explicit TrialBeta(argtype args = argtype());
  explicit TrialBeta(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialBeta>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialBeta>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialBeta(std::istream& istr);
  virtual ~TrialBeta() {}
  //@}
};

inline std::shared_ptr<TrialBeta> MakeTrialBeta(argtype args = argtype()) {
  return std::make_shared<TrialBeta>(args); }

}  // namespace feasst

#endif  // FEASST_BETA_EXPANDED_TRIAL_BETA_H_

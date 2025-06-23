#ifndef FEASST_MONTE_CARLO_TRIAL_ADD_H_
#define FEASST_MONTE_CARLO_TRIAL_ADD_H_

#include <memory>
#include "monte_carlo/include/trial.h"

namespace feasst {

/**
 * Attempt to add a particle to intialize a Configuration.
 * This trial does not obey detailed balance without a TrialRemove of the same
 * weight.
 * Instead, use TrialTransfer to ensure detailed balance.
 */
class TrialAdd : public Trial {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - TrialSelectParticle arguments.
    - TrialStage arguments.
   */
  explicit TrialAdd(argtype args = argtype());
  explicit TrialAdd(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAdd>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAdd>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAdd(std::istream& istr);
  virtual ~TrialAdd() {}
  //@}
};

inline std::shared_ptr<TrialAdd> MakeTrialAdd(argtype args = argtype()) {
  return std::make_shared<TrialAdd>(args); }

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_ADD_H_


#ifndef FEASST_CHAIN_TRIAL_SWAP_POSITION_H_
#define FEASST_CHAIN_TRIAL_SWAP_POSITION_H_

#include <map>
#include <string>
#include <memory>
#include "monte_carlo/include/trial_move.h"

namespace feasst {
/**
 *
 */
class TrialSwapPosition : public TrialMove {
 public:
  //@{
  /** @name Arguments
    args:
    - TrialSelectParticle arguments.
    - TrialStage arguments.
    - Tunable arguments.
   */
  explicit TrialSwapPosition(argtype args = argtype());
  explicit TrialSwapPosition(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialSwapPosition>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialSwapPosition>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialSwapPosition(std::istream& istr);
  virtual ~TrialSwapPosition() {}
};

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SWAP_POSITION_H_

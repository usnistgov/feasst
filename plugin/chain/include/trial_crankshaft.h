#ifndef FEASST_CHAIN_TRIAL_CRANKSHAFT_H_
#define FEASST_CHAIN_TRIAL_CRANKSHAFT_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

/// Rigidly rotate a sub section of a chain.
class TrialCrankshaft : public TrialMove {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - TrialStage arguments.
    - SelectSegment arguments.
    - Tunable arguments.
   */
  explicit TrialCrankshaft(argtype args = argtype());
  explicit TrialCrankshaft(argtype * args);
  //@}
  /** @name Public Functions
   */
  //@{
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialCrankshaft>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialCrankshaft>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialCrankshaft(std::istream& istr);
  virtual ~TrialCrankshaft() {}
  //@}
};

inline std::shared_ptr<TrialCrankshaft> MakeTrialCrankshaft(argtype args = argtype()) {
  return std::make_shared<TrialCrankshaft>(args); }

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_CRANKSHAFT_H_

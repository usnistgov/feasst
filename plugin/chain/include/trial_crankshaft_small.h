#ifndef FEASST_CHAIN_TRIAL_CRANKSHAFT_SMALL_H_
#define FEASST_CHAIN_TRIAL_CRANKSHAFT_SMALL_H_

#include <memory>
#include "monte_carlo/include/trial_move.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Rigidly rotate a sub section of a chain.
  This Trial is optimized for smaller molecules.
 */
class TrialCrankshaftSmall : public TrialMove {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - TrialStage arguments.
    - SelectCrankshaftSmall arguments.
    - Tunable arguments.
   */
  explicit TrialCrankshaftSmall(argtype args = argtype());
  explicit TrialCrankshaftSmall(argtype * args);
  //@}
  /** @name Public Functions
   */
  //@{
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialCrankshaftSmall>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialCrankshaftSmall>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialCrankshaftSmall(std::istream& istr);
  virtual ~TrialCrankshaftSmall() {}
  //@}
};

inline std::shared_ptr<TrialCrankshaftSmall> MakeTrialCrankshaftSmall(argtype args = argtype()) {
  return std::make_shared<TrialCrankshaftSmall>(args); }

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_CRANKSHAFT_SMALL_H_

#ifndef FEASST_CHAIN_TRIAL_CRANKSHAFT_SMALL_H_
#define FEASST_CHAIN_TRIAL_CRANKSHAFT_SMALL_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

/// Rigidly rotate a sub section of a chain.
class TrialCrankshaftSmall : public TrialMove {
 public:
  TrialCrankshaftSmall(argtype args = argtype());
  TrialCrankshaftSmall(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialCrankshaftSmall>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialCrankshaftSmall>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialCrankshaftSmall(std::istream& istr);
  virtual ~TrialCrankshaftSmall() {}
};

inline std::shared_ptr<TrialCrankshaftSmall> MakeTrialCrankshaftSmall(argtype args = argtype()) {
  return std::make_shared<TrialCrankshaftSmall>(args); }

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_CRANKSHAFT_SMALL_H_

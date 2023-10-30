#ifndef FEASST_CHAIN_TRIAL_PIVOT_H_
#define FEASST_CHAIN_TRIAL_PIVOT_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

/// Rigidly pivot an end segment of a chain to a random orientation.
class TrialPivot : public TrialMove {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - TrialStage arguments.
    - SelectEndSegment arguments.
    - Tunable arguments.
   */
  explicit TrialPivot(argtype args = argtype());
  explicit TrialPivot(argtype * args);
  //@}
  /** @name Public Functions
   */
  //@{
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialPivot>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialPivot>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialPivot(std::istream& istr);
  virtual ~TrialPivot() {}
  //@}
};

inline std::shared_ptr<TrialPivot> MakeTrialPivot(argtype args = argtype()) {
  return std::make_shared<TrialPivot>(args); }

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_PIVOT_H_

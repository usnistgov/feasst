
#ifndef FEASST_CHAIN_TRIAL_PIVOT_H_
#define FEASST_CHAIN_TRIAL_PIVOT_H_

#include <memory>
#include "monte_carlo/include/trial_move.h"
#include "chain/include/trial_select_end_segment.h"
#include "chain/include/perturb_pivot.h"

namespace feasst {

class TrialPivot : public TrialMove {
 public:
  TrialPivot(
    /// These arguments are sent to both PerturbPivot and TrialStage.
    const argtype& args = argtype())
    : TrialMove(
      std::make_shared<TrialSelectEndSegment>(args),
      std::make_shared<PerturbPivot>(args),
      args
    ) {
    class_name_ = "TrialPivot";
  }
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialPivot(std::istream& istr);
  virtual ~TrialPivot() {}
};

inline std::shared_ptr<TrialPivot> MakeTrialPivot(
    const argtype &args = argtype()) {
  return std::make_shared<TrialPivot>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_PIVOT_H_

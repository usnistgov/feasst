
#ifndef FEASST_CHAIN_TRIAL_H_
#define FEASST_CHAIN_TRIAL_H_

#include <string>
#include <memory>
#include "monte_carlo/include/trial_move.h"
#include "chain/include/trial_select_reptate.h"
#include "chain/include/perturb_reptate.h"

namespace feasst {

class TrialReptate : public TrialMove {
 public:
  TrialReptate(
    /// These arguments are sent to both PerturbReptate and TrialStage.
    const argtype& args = argtype())
    : TrialMove(
      std::make_shared<TrialSelectReptate>(),
      std::make_shared<PerturbReptate>(args),
      args
    ) {
    class_name_ = "TrialReptate";
  }
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialReptate(std::istream& istr);
  virtual ~TrialReptate() {}
};

inline std::shared_ptr<TrialReptate> MakeTrialReptate(
    const argtype &args = argtype()) {
  return std::make_shared<TrialReptate>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_H_

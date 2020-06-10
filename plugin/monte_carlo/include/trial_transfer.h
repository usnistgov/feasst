
#ifndef FEASST_MONTE_CARLO_TRIAL_TRANSFER_H_
#define FEASST_MONTE_CARLO_TRIAL_TRANSFER_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Attempt TrialAdd or TrialRemove with equal probability.
class TrialTransfer : public TrialFactory {
 public:
  explicit TrialTransfer(const argtype& args = argtype());
};

inline std::shared_ptr<TrialTransfer> MakeTrialTransfer(
    const argtype &args = argtype()) {
  return std::make_shared<TrialTransfer>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_TRANSFER_H_

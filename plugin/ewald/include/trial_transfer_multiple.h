
#ifndef FEASST_MONTE_CARLO_TRIAL_TRANSFER_MULTIPLE_H_
#define FEASST_MONTE_CARLO_TRIAL_TRANSFER_MULTIPLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Attempt TrialAddMultiple or TrialRemoveMultiple with equal probability.
class TrialTransferMultiple : public TrialFactory {
 public:
  explicit TrialTransferMultiple(const argtype& args = argtype());
};

inline std::shared_ptr<TrialTransferMultiple> MakeTrialTransferMultiple(
    const argtype &args = argtype()) {
  return std::make_shared<TrialTransferMultiple>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_TRANSFER_MULTIPLE_H_

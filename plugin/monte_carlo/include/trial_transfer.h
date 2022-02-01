#ifndef FEASST_MONTE_CARLO_TRIAL_TRANSFER_H_
#define FEASST_MONTE_CARLO_TRIAL_TRANSFER_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Add both a TrialAdd and TrialRemove with the same arguments.
class TrialTransfer : public TrialFactoryNamed {
 public:
  TrialTransfer(argtype args = argtype());
  TrialTransfer(argtype * args);
  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialTransfer>(args); }
  virtual ~TrialTransfer() {}
};

inline std::shared_ptr<TrialTransfer> MakeTrialTransfer(argtype args = argtype()) {
  return std::make_shared<TrialTransfer>(args); }

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_TRANSFER_H_

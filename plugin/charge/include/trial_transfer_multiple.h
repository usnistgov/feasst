
#ifndef FEASST_CHARGE_TRIAL_TRANSFER_MULTIPLE_H_
#define FEASST_CHARGE_TRIAL_TRANSFER_MULTIPLE_H_

#include <memory>
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/// Attempt TrialAddMultiple or TrialRemoveMultiple and split the trial weights
/// equally.
class TrialTransferMultiple : public TrialFactoryNamed {
 public:
  //@{
  /** @name Arguments
    - particle_types: comma-separated list of particle names to transfer.
    - TrialStage arguments.
   */
  explicit TrialTransferMultiple(argtype args = argtype());
  explicit TrialTransferMultiple(argtype * args);
  //@}
  /** @name Public Functions
   */
  //@{
  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialTransferMultiple>(args); }
  virtual ~TrialTransferMultiple() {}
  //@}
};

inline std::shared_ptr<TrialTransferMultiple> MakeTrialTransferMultiple(argtype args = argtype()) {
  return std::make_shared<TrialTransferMultiple>(args); }

}  // namespace feasst

#endif  // FEASST_CHARGE_TRIAL_TRANSFER_MULTIPLE_H_


#ifndef FEASST_CLUSTER_TRIAL_TRANSFER_AVB_H_
#define FEASST_CLUSTER_TRIAL_TRANSFER_AVB_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Attempt TrialAddAVB or TrialRemoveAVB with equal probability.
class TrialTransferAVB : public TrialFactoryNamed {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - SelectParticleAVB arguments.
   */
  explicit TrialTransferAVB(argtype args = argtype());
  explicit TrialTransferAVB(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialTransferAVB>(args); }
  virtual ~TrialTransferAVB() {}
  //@}
};

inline std::shared_ptr<TrialTransferAVB> MakeTrialTransferAVB(argtype args = argtype()) {
  return std::make_shared<TrialTransferAVB>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_TRANSFER_AVB_H_

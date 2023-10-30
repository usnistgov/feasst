
#ifndef FEASST_CLUSTER_TRIAL_TRANSFER_AVB_DIVALENT_H_
#define FEASST_CLUSTER_TRIAL_TRANSFER_AVB_DIVALENT_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Attempt TrialAddAVBDivalent or TrialRemoveAVBDivalent with equal probability
class TrialTransferAVBDivalent : public TrialFactoryNamed {
 public:
  //@{
  /** @name Arguments
    - particle_type_a: type of second added particle in AV of first.
    - site_index_a: index of site in type a that defines AV (default: 0).
    - particle_type_b: type of third added particle in AV of first.
    - site_index_b: index of site in type b that defines AV (default: 0).
    - SelectParticleAVBDivalent arguments.
    - Trial arguments.
   */
  explicit TrialTransferAVBDivalent(argtype args = argtype());
  explicit TrialTransferAVBDivalent(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialTransferAVBDivalent>(args); }
  virtual ~TrialTransferAVBDivalent() {}
  //@}
};

inline std::shared_ptr<TrialTransferAVBDivalent> MakeTrialTransferAVBDivalent(argtype args = argtype()) {
  return std::make_shared<TrialTransferAVBDivalent>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_TRANSFER_AVB_DIVALENT_H_

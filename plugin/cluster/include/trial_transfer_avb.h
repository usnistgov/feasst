
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
   */

  /**
    args:
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

/// Attempt to add a particle with AVB as described in ComputeAddAVB.
class TrialAddAVB : public Trial {
 public:
  explicit TrialAddAVB(argtype args = argtype());
  explicit TrialAddAVB(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAddAVB>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAddAVB>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAddAVB(std::istream& istr);
  virtual ~TrialAddAVB() {}
};

inline std::shared_ptr<TrialAddAVB> MakeTrialAddAVB(argtype args = argtype()) {
  return std::make_shared<TrialAddAVB>(args); }

/// Attempt to remove a particle with AVB as described in ComputeRemoveAVB.
class TrialRemoveAVB : public Trial {
 public:
  explicit TrialRemoveAVB(argtype args = argtype());
  explicit TrialRemoveAVB(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialRemoveAVB>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialRemoveAVB>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialRemoveAVB(std::istream& istr);
  virtual ~TrialRemoveAVB() {}
};

inline std::shared_ptr<TrialRemoveAVB> MakeTrialRemoveAVB(argtype args = argtype()) {
  return std::make_shared<TrialRemoveAVB>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_TRANSFER_AVB_H_

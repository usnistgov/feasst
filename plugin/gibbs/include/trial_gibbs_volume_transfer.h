
#ifndef FEASST_GIBBS_TRIAL_GIBBS_VOLUME_TRANSFER_H_
#define FEASST_GIBBS_TRIAL_GIBBS_VOLUME_TRANSFER_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Attempt TrialGibbsVolumeTransferOneWay with equal probability in either direction.
  See https://doi.org/10.1080/00268978700101491.
 */
class TrialGibbsVolumeTransfer : public TrialFactoryNamed {
 public:
  //@{
  /** @name Arguments
    - configuration_index0: index of one of the configurations (default: 0).
    - configuration_index1: index of the other configuration (default: 1).
    - Tunable arguments.
   */
  explicit TrialGibbsVolumeTransfer(argtype args = argtype());
  explicit TrialGibbsVolumeTransfer(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialGibbsVolumeTransfer>(args); }
  virtual ~TrialGibbsVolumeTransfer() {}
  //@}
};

inline std::shared_ptr<TrialGibbsVolumeTransfer> MakeTrialGibbsVolumeTransfer(argtype args = argtype()) {
  return std::make_shared<TrialGibbsVolumeTransfer>(args); }

/// Attempt to transfer volume from one configuration to another.
class TrialGibbsVolumeTransferOneWay : public Trial {
 public:
  /**
    args:
    - to_configuration_index: index of configuration to send the particle.
    - configuration_index: from TrialSelect, configuration which donates a
      particle (default: 0).
   */
  explicit TrialGibbsVolumeTransferOneWay(argtype args = argtype());
  explicit TrialGibbsVolumeTransferOneWay(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialGibbsVolumeTransferOneWay>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialGibbsVolumeTransferOneWay>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialGibbsVolumeTransferOneWay(std::istream& istr);
  virtual ~TrialGibbsVolumeTransferOneWay() {}
};

inline std::shared_ptr<TrialGibbsVolumeTransferOneWay> MakeTrialGibbsVolumeTransferOneWay(argtype args = argtype()) {
  return std::make_shared<TrialGibbsVolumeTransferOneWay>(args); }

}  // namespace feasst

#endif  // FEASST_GIBBS_TRIAL_GIBBS_VOLUME_TRANSFER_H_

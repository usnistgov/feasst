
#ifndef FEASST_GIBBS_TRIAL_GIBBS_VOLUME_TRANSFER_H_
#define FEASST_GIBBS_TRIAL_GIBBS_VOLUME_TRANSFER_H_

#include <memory>
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Attempt to transfer a uniformly random amount of volume between two
  configurations.
  See https://doi.org/10.1080/00268978700101491.
 */
class TrialGibbsVolumeTransfer : public Trial {
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

  // Same as Trial but also check reference potential.
  void precompute(Criteria * criteria, System * system) override;

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialGibbsVolumeTransfer>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialGibbsVolumeTransfer>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialGibbsVolumeTransfer(std::istream& istr);
  virtual ~TrialGibbsVolumeTransfer() {}
  //@}
};

}  // namespace feasst

#endif  // FEASST_GIBBS_TRIAL_GIBBS_VOLUME_TRANSFER_H_

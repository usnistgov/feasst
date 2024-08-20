
#ifndef FEASST_GIBBS_TRIAL_GIBBS_PARTICLE_TRANSFER_H_
#define FEASST_GIBBS_TRIAL_GIBBS_PARTICLE_TRANSFER_H_

#include <memory>
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Attempt TrialGibbsParticleTransferOneWay with equal probability in either direction.
  See https://doi.org/10.1080/00268978700101491.
 */
class TrialGibbsParticleTransfer : public TrialFactoryNamed {
 public:
  //@{
  /** @name Arguments
    - configuration_index0: index of one of the configurations (default: 0).
    - configuration_index1: index of the other configuration (default: 1).
    - TrialSelectParticle arguments.
   */
  explicit TrialGibbsParticleTransfer(argtype args = argtype());
  explicit TrialGibbsParticleTransfer(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialGibbsParticleTransfer>(args); }
  virtual ~TrialGibbsParticleTransfer() {}
  //@}
};

inline std::shared_ptr<TrialGibbsParticleTransfer> MakeTrialGibbsParticleTransfer(argtype args = argtype()) {
  return std::make_shared<TrialGibbsParticleTransfer>(args); }

/// Attempt to transfer a particle from one configuration to another.
class TrialGibbsParticleTransferOneWay : public Trial {
 public:
  /**
    args:
    - to_configuration_index: index of configuration to send the particle.
    - configuration_index: from TrialSelect, configuration which donates a
      particle (default: 0).
   */
  explicit TrialGibbsParticleTransferOneWay(argtype args = argtype());
  explicit TrialGibbsParticleTransferOneWay(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialGibbsParticleTransferOneWay>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialGibbsParticleTransferOneWay>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialGibbsParticleTransferOneWay(std::istream& istr);
  virtual ~TrialGibbsParticleTransferOneWay() {}
};

inline std::shared_ptr<TrialGibbsParticleTransferOneWay> MakeTrialGibbsParticleTransferOneWay(argtype args = argtype()) {
  return std::make_shared<TrialGibbsParticleTransferOneWay>(args); }

}  // namespace feasst

#endif  // FEASST_GIBBS_TRIAL_GIBBS_PARTICLE_TRANSFER_H_


#ifndef FEASST_GIBBS_COMPUTE_GIBBS_PARTICLE_TRANSFER_H_
#define FEASST_GIBBS_COMPUTE_GIBBS_PARTICLE_TRANSFER_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

/**
  Attempt to transfer a particle from one configuration to another,
  and vice versa.
 */
class ComputeGibbsParticleTransfer : public TrialCompute {
 public:
  ComputeGibbsParticleTransfer();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeGibbsParticleTransfer(std::istream& istr);
  virtual ~ComputeGibbsParticleTransfer() {}

 protected:
  void serialize_compute_gibbs_particle_transfer_(std::ostream& ostr) const;
};

inline std::shared_ptr<ComputeGibbsParticleTransfer> MakeComputeGibbsParticleTransfer() {
  return std::make_shared<ComputeGibbsParticleTransfer>();
}

}  // namespace feasst

#endif  // FEASST_GIBBS_COMPUTE_GIBBS_PARTICLE_TRANSFER_H_

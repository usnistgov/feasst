
#ifndef FEASST_GIBBS_COMPUTE_GIBBS_VOLUME_TRANSFER_H_
#define FEASST_GIBBS_COMPUTE_GIBBS_VOLUME_TRANSFER_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

/**
  Attempt to transfer volome from one configuration to another,
  and vice versa.
 */
class ComputeGibbsVolumeTransfer : public TrialCompute {
 public:
  ComputeGibbsVolumeTransfer();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeGibbsVolumeTransfer(std::istream& istr);
  virtual ~ComputeGibbsVolumeTransfer() {}

 protected:
  void serialize_compute_gibbs_volume_transfer_(std::ostream& ostr) const;
};

inline std::shared_ptr<ComputeGibbsVolumeTransfer> MakeComputeGibbsVolumeTransfer() {
  return std::make_shared<ComputeGibbsVolumeTransfer>();
}

}  // namespace feasst

#endif  // FEASST_GIBBS_COMPUTE_GIBBS_VOLUME_TRANSFER_H_

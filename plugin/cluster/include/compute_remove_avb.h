
#ifndef FEASST_MONTE_CARLO_COMPUTE_REMOVE_AVB_H_
#define FEASST_MONTE_CARLO_COMPUTE_REMOVE_AVB_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

/**
  Attempt to remove a particle of type "a" with a bias for site_index_a
  to be in the aggrevation volume (AV) of site_index_t ("t" for target) of
  particle of type "t".

  See ComputeAddAVB for derivation of the acceptance probability that is
  the reverse of this Trial.

  There are modifications to make for this reverse move considering that the
  acceptance probability is computed before the removal takes place.

  The number of sites to select in the AV already contains the site added from
  the old state,

  \f$ N^{s,AV}_a + 1 \rightarrow N^{s,AV}_a \f$.
 */
class ComputeRemoveAVB : public TrialCompute {
 public:
  ComputeRemoveAVB();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeRemoveAVB(std::istream& istr);
  virtual ~ComputeRemoveAVB() {}

 protected:
  void serialize_compute_remove_avb_(std::ostream& ostr) const;
};

inline std::shared_ptr<ComputeRemoveAVB> MakeComputeRemoveAVB() {
  return std::make_shared<ComputeRemoveAVB>();
}
}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_COMPUTE_REMOVE_AVB_H_

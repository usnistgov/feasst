
#ifndef FEASST_CLUSTER_COMPUTE_REMOVE_AVB_DIVALENT_H_
#define FEASST_CLUSTER_COMPUTE_REMOVE_AVB_DIVALENT_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

/**
  Attempt to remove a particle of type "t" with site_index_t anywhere in Domain.
  Then, remove a second particle of type "a" with site_index_a in the AV of
  site_index_t.
  Finally, remove a third particle of type "a" with site_index_a in the AV of
  site_index_t.

  See ComputeAddAVBDivalent for derivation of the acceptance probability that
  is the reverse of this Trial.

  There are modifications to make for this reverse move considering that the
  acceptance probability is computed before the removal takes place.

  The number of sites to select in the AV already contains the site added from
  the old state,

  \f$ N^{s,AV}_a + [1,2] \rightarrow N^{s,AV}_a \f$.

  The number of particles of type "t" already contains the first particle added
  from the old state.

  \f$ N_t + 1 \rightarrow N_t \f$
 */
class ComputeRemoveAVBDivalent : public TrialCompute {
 public:
  explicit ComputeRemoveAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria);

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeRemoveAVBDivalent(std::istream& istr);
  virtual ~ComputeRemoveAVBDivalent() {}

 protected:
  void serialize_compute_remove_avb_divalent_(std::ostream& ostr) const;

 private:
  std::shared_ptr<NeighborCriteria> neighbor_criteria_;

  // temporary
  Select neighbors_;
};

inline std::shared_ptr<ComputeRemoveAVBDivalent> MakeComputeRemoveAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria) {
  return std::make_shared<ComputeRemoveAVBDivalent>(neighbor_criteria);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_COMPUTE_REMOVE_AVB_DIVALENT_H_

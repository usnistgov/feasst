
#ifndef FEASST_CLUSTER_COMPUTE_ADD_AVB_DIVALENT_H_
#define FEASST_CLUSTER_COMPUTE_ADD_AVB_DIVALENT_H_

#include <memory>
#include <vector>
#include "configuration/include/select.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class ComputeAddAVBDivalent : public TrialCompute {
 public:
  /**
    args:
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
   */
  explicit ComputeAddAVBDivalent(argtype args = argtype());

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeAddAVBDivalent(std::istream& istr);
  virtual ~ComputeAddAVBDivalent() {}

 protected:
  void serialize_compute_add_avb_divalent_(std::ostream& ostr) const;

 private:
  int neighbor_;

  // temporary
  Select neighbors_;
};

inline std::shared_ptr<ComputeAddAVBDivalent> MakeComputeAddAVBDivalent(
    argtype args = argtype()) {
  return std::make_shared<ComputeAddAVBDivalent>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_COMPUTE_ADD_AVB_DIVALENT_H_

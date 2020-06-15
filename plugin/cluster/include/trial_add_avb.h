
#ifndef FEASST_CLUSTER_TRIAL_ADD_AVB_H_
#define FEASST_CLUSTER_TRIAL_ADD_AVB_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "system/include/neighbor_criteria.h"

namespace feasst {

/// Attempt to add a particle with AVB as described in ComputeAddAVB.
class TrialAddAVB : public Trial {
 public:
  TrialAddAVB(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialAddAVB(std::istream& istr);
  virtual ~TrialAddAVB() {}

 protected:
  void serialize_trial_add_avb_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialAddAVB> MakeTrialAddAVB(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype()) {
  return std::make_shared<TrialAddAVB>(neighbor_criteria, args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_ADD_AVB_H_

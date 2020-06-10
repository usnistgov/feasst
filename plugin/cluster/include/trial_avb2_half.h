
#ifndef FEASST_CLUSTER_TRIAL_AVB2_HALF_H_
#define FEASST_CLUSTER_TRIAL_AVB2_HALF_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "system/include/neighbor_criteria.h"

namespace feasst {

class TrialAVB2Half : public Trial {
 public:
  TrialAVB2Half(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialAVB2Half(std::istream& istr);
  virtual ~TrialAVB2Half() {}

 protected:
  void serialize_trial_avb2_half_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialAVB2Half> MakeTrialAVB2Half(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype()) {
  return std::make_shared<TrialAVB2Half>(neighbor_criteria, args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_AVB2_HALF_H_

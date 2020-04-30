
#ifndef FEASST_CLUSTER_TRIAL_AVB4_H_
#define FEASST_CLUSTER_TRIAL_AVB4_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "system/include/neighbor_criteria.h"

namespace feasst {

class TrialAVB4 : public Trial {
 public:
  TrialAVB4(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialAVB4(std::istream& istr);
  virtual ~TrialAVB4() {}

 protected:
  void serialize_trial_avb4_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialAVB4> MakeTrialAVB4(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype()) {
  return std::make_shared<TrialAVB4>(neighbor_criteria, args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_AVB4_H_

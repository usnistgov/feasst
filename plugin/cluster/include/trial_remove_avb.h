
#ifndef FEASST_CLUSTER_TRIAL_REMOVE_AVB_H_
#define FEASST_CLUSTER_TRIAL_REMOVE_AVB_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Attempt to remove a particle with AVB as described in ComputeRemoveAVB.
class TrialRemoveAVB : public Trial {
 public:
  explicit TrialRemoveAVB(argtype args = argtype());
  explicit TrialRemoveAVB(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialRemoveAVB>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialRemoveAVB>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialRemoveAVB(std::istream& istr);
  virtual ~TrialRemoveAVB() {}
};

inline std::shared_ptr<TrialRemoveAVB> MakeTrialRemoveAVB(argtype args = argtype()) {
  return std::make_shared<TrialRemoveAVB>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_REMOVE_AVB_H_

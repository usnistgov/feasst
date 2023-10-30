
#ifndef FEASST_CLUSTER_TRIAL_ADD_AVB_H_
#define FEASST_CLUSTER_TRIAL_ADD_AVB_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Attempt to add a particle with AVB as described in ComputeAddAVB.
class TrialAddAVB : public Trial {
 public:
  explicit TrialAddAVB(argtype args = argtype());
  explicit TrialAddAVB(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAddAVB>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAddAVB>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAddAVB(std::istream& istr);
  virtual ~TrialAddAVB() {}
};

inline std::shared_ptr<TrialAddAVB> MakeTrialAddAVB(argtype args = argtype()) {
  return std::make_shared<TrialAddAVB>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_ADD_AVB_H_

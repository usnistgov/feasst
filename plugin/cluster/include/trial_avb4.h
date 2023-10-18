
#ifndef FEASST_CLUSTER_TRIAL_AVB4_H_
#define FEASST_CLUSTER_TRIAL_AVB4_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/**
  Only implemented for single site particles.
  See ComputeAVB4.
 */
class TrialAVB4 : public Trial {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
    - Trial arguments.
    - TrialStage arguments.
    - SelectParticleAVB arguments.
   */
  explicit TrialAVB4(argtype args = argtype());
  explicit TrialAVB4(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAVB4>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAVB4>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAVB4(std::istream& istr);
  virtual ~TrialAVB4() {}
  //@}
};

inline std::shared_ptr<TrialAVB4> MakeTrialAVB4(argtype args = argtype()) {
  return std::make_shared<TrialAVB4>(args); }

void gen_avb4_args_(argtype * args);

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_AVB4_H_

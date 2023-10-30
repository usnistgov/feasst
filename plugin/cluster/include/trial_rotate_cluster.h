
#ifndef FEASST_CLUSTER_TRIAL_ROTATE_CLUSTER_H_
#define FEASST_CLUSTER_TRIAL_ROTATE_CLUSTER_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Rigid rotation of a cluster of particles.
 */
class TrialRotateCluster : public Trial {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - SelectCluster arguments.
    - Tunable arguments.
   */
  explicit TrialRotateCluster(argtype args = argtype());
  explicit TrialRotateCluster(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialRotateCluster>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialRotateCluster>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialRotateCluster(std::istream& istr);
  virtual ~TrialRotateCluster() {}
  //@}
};

inline std::shared_ptr<TrialRotateCluster> MakeTrialRotateCluster(argtype args = argtype()) {
  return std::make_shared<TrialRotateCluster>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_ROTATE_CLUSTER_H_


#ifndef FEASST_CLUSTER_TRIAL_RIGID_CLUSTER_H_
#define FEASST_CLUSTER_TRIAL_RIGID_CLUSTER_H_

#include <memory>
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Attempt TrialTranslateCluster and TrialRotateCluster with equal probability.
 */
class TrialRigidCluster : public TrialFactoryNamed {
 public:
  //@{
  /** @name Arguments
    - rotate_param: initial value of the tunable parameter (default: 25).
    - translate_param: initial value of the tunable parameter (default: 0.1).
    - Trial arguments.
    - SelectCluster arguments.
    - Tunable arguments.
   */
  explicit TrialRigidCluster(argtype args = argtype());
  explicit TrialRigidCluster(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialRigidCluster>(args); }
  virtual ~TrialRigidCluster() {}
  //@}
};

inline std::shared_ptr<TrialRigidCluster> MakeTrialRigidCluster(argtype args = argtype()) {
  return std::make_shared<TrialRigidCluster>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_RIGID_CLUSTER_H_

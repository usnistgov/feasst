
#ifndef FEASST_CLUSTER_TRIAL_RIGID_CLUSTER_H_
#define FEASST_CLUSTER_TRIAL_RIGID_CLUSTER_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Rigid translation of a cluster of particles.
 */
std::shared_ptr<Trial> MakeTrialTranslateCluster(
  const argtype &args = argtype());

/**
  Rigid translation of a cluster of particles.
 */
std::shared_ptr<Trial> MakeTrialRotateCluster(
  const argtype &args = argtype());

/**
  Attempt TrialTranslateCluster and TrialRotateCluster with equal probability.

  args:
  - rotate_param: initial value of the tunable parameter (default: 25).
  - translate_param: initial value of the tunable parameter (default: 0.1).
 */
std::shared_ptr<TrialFactory> MakeTrialRigidCluster(
  const argtype &args = argtype());

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_RIGID_CLUSTER_H_

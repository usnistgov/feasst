
#ifndef FEASST_CLUSTER_UTILS_CLUSTER_H_
#define FEASST_CLUSTER_UTILS_CLUSTER_H_

#include "utils/include/arguments.h"
#include "monte_carlo/include/monte_carlo.h"
#include "cluster/include/trial_rigid_cluster.h"

namespace feasst {

/// Initialize a cluster translation and rotation trial simultaneously with the
/// same arguments.
inline void add_rigid_cluster_trials(MonteCarlo * mc,
    std::shared_ptr<ClusterCriteria> cluster_criteria,
    const argtype& args = argtype()) {
  mc->add(MakeTrialTranslateCluster(cluster_criteria, args));
  mc->add(MakeTrialRotateCluster(cluster_criteria, args));
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_UTILS_CLUSTER_H_

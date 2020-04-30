
#ifndef FEASST_CLUSTER_UTILS_CLUSTER_H_
#define FEASST_CLUSTER_UTILS_CLUSTER_H_

#include "utils/include/arguments.h"
#include "monte_carlo/include/monte_carlo.h"
#include "cluster/include/trial_rigid_cluster.h"
#include "cluster/include/trial_add_avb.h"
#include "cluster/include/trial_remove_avb.h"
#include "cluster/include/trial_avb2.h"

namespace feasst {

/// Initialize a cluster translation and rotation trial simultaneously with the
/// same arguments.
inline void add_rigid_cluster_trials(MonteCarlo * mc,
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype()) {
  mc->add(MakeTrialTranslateCluster(neighbor_criteria, args));
  mc->add(MakeTrialRotateCluster(neighbor_criteria, args));
}

/// Add an AVB insertion/deletion trial to mc.
inline void add_avb_transfer_trials(MonteCarlo * mc,
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype()) {
  mc->add(MakeTrialAddAVB(neighbor_criteria, args));
  mc->add(MakeTrialRemoveAVB(neighbor_criteria, args));
}

/// Add AVB2 trials to mc.
inline void add_avb2_trials(MonteCarlo *mc,
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype()) {
  argtype args1(args);
  args1.insert({{"out_to_in", "true"}});
  mc->add(MakeTrialAVB2(neighbor_criteria, args1));
  argtype args2(args);
  args2.insert({{"out_to_in", "false"}});
  mc->add(MakeTrialAVB2(neighbor_criteria, args2));
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_UTILS_CLUSTER_H_

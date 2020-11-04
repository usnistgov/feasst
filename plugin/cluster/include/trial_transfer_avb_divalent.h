
#ifndef FEASST_CLUSTER_TRIAL_ADD_AVB_DIVALENT_H_
#define FEASST_CLUSTER_TRIAL_ADD_AVB_DIVALENT_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"
#include "system/include/neighbor_criteria.h"

namespace feasst {

/**
  Attempt to add a particle with AVB as described in ComputeAddAVBDivalent.

  args:
  - particle_type_a: type of second added particle in AV of first.
  - site_index_a: index of site in type a that defines AV (default: 0).
  - particle_type_b: type of third added particle in AV of first.
  - site_index_b: index of site in type b that defines AV (default: 0).
 */
std::shared_ptr<Trial> MakeTrialAddAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype());

/**
  Attempt to add a particle with AVB as described in ComputeRemoveAVBDivalent.

  args:
  - particle_type_a: type of second added particle in AV of first.
  - site_index_a: index of site in type a that defines AV (default: 0).
  - particle_type_b: type of third added particle in AV of first.
  - site_index_b: index of site in type b that defines AV (default: 0).
 */
std::shared_ptr<Trial> MakeTrialRemoveAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype());

/// Attempt TrialAddAVBDivalent or TrialRemoveAVBDivalent with equal probability
std::shared_ptr<TrialFactory> MakeTrialTransferAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype());

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_ADD_AVB_DIVALENT_H_

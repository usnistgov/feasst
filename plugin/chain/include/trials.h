#ifndef FEASST_CHAIN_TRIALS_H_
#define FEASST_CHAIN_TRIALS_H_

#include <memory>
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Rigidly pivot an end segment of a chain to a random orientation.
std::shared_ptr<Trial> MakeTrialPivot(const argtype &args = argtype());

/// Rigidly rotate a sub section of a chain.
std::shared_ptr<Trial> MakeTrialCrankshaft(const argtype &args = argtype());

/**
  Reptate a linear chain by taking one end and adding it to the other end.
  For heteropolymers, this perturbation changes the composition of all
  particles with the same type.
  Thus, individual heteropolymers should be added as unique particle types.
 */
std::shared_ptr<Trial> MakeTrialReptate(const argtype &args = argtype());

/**
  Swap the types of two sites in a particle.
  args:
  - site_type1: type of site to swap.
  - site_type2: type of other site to swap.
 */
std::shared_ptr<Trial> MakeTrialSwapSites(const argtype &args = argtype());

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIALS_H_


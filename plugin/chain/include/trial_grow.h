#ifndef FEASST_CHAIN_TRIAL_GROW_H_
#define FEASST_CHAIN_TRIAL_GROW_H_

#include <memory>
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Manually describe partial or full-particle growth using configurational bias.
  The input is a vector of argtype, where each argtype represents a TrialStage
  for growing one (or rarely, multiple) sites.

  The following options may only be used in the first argtype.
  - particle_type: type of particle in Configuration (always required).
  - site: site index in particle_type to transfer/regrow (always required).
  - weight: weight of selection of this trial (default: see Trial).
  - transfer: if true, create add and remove trial with equal weight
    (default: false).
  - regrow: if true, place anywhere in the domain (default: false).
  - transfer_avb: if true, same as transfer but with AVB Add/Remove for the
    first stage (default: false).
  - regrow_avb2: if true, regrow using AVB2 in the first stage (default: false).
  - regrow_avb4: if true, regrow using AVB4 in the first stage (default: false).

  The following options may be used in any argtype, including the first for
  particle regrowth.
  - num_steps: number of Rosenbluth steps (see TrialStage).
  - reference_index: reference potential for DCCB (see TrialStage).
  - bond: if used, add TrialSelectBond and PerturbDistance.
    Requires arguments described in TrialSelectBond.
  - angle: if used, adds TrialSelectAngle and PerturbDistanceAngle.
    Requires arguments described in TrialSelectAngle and TrialSelectBond.
  - branch: if used, adds SelectBranch and PerturbBranch.
    Requires arguments described in SelectBranch, Angle and Bond.

  Note that only one of bond, angle or branch may be true for a given stage.
 */
std::shared_ptr<TrialFactory> MakeTrialGrow(
    const std::vector<argtype>& args,
    /// Optionally, the following arguments for num_steps and reference_index
    /// are applied to every growth.
    /// Any option applied by the above args overwrites this option.
    const argtype& default_args = argtype());

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_GROW_H_

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
  - translate: if true (default: false), translate site (which is required arg
    for TrialSelectParticle).
    In addition, must have number of stages equal to number of sites.

  The following options may be used in any argtype.
  If used in the first, then its a partial regrowth move.
  - bond: if true, add TrialSelectBond and PerturbDistance.
    Requires arguments described in TrialSelectBond.
  - angle: if true, adds TrialSelectAngle and PerturbDistanceAngle.
    Requires arguments described in TrialSelectAngle and TrialSelectBond.
  - dihedral: if true, adds TrialSelectDihedral and PerturbDihedral.
    Requires arguments described in SelectDihedral, Angle and Bond.
  - branch: if true, adds SelectBranch and PerturbBranch.
    Requires arguments described in SelectBranch, Angle and Bond.
  - reptate: if true, add TrialSelectBond and PerturbToAnchor.
    Requires arguments described in TrialSelectBond.
  - TrialStage arguments: num_steps, reference_index, new_only, etc.

  Note that only one of bond, angle or branch may be true for a given stage.

  Derivation of the acceptance criteria follows the procedure described in
  TrialComputeMove and first described in Lecture 9 of Prof. David Kofke's
  Molecular Simulation course CE 530

  http://www.eng.buffalo.edu/~kofke/ce530/Lectures/lectures.html

  For bonds, angles and dihedrals, the relevant forward and reverse
  probabilities are as follows:
\rst
+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Select particle of type t            |:math:`1/N_t`                           |
|                                     |                                        |
|[Select particle of type t]          |:math:`[1/N_t]`                         |
+-------------------------------------+----------------------------------------+
|Select bond from given site indices  |:math:`1`                               |
|                                     |                                        |
|[Select bond from given site indices]|:math:`1`                               |
+-------------------------------------+----------------------------------------+
|Select new length/angle :math:`b_n`  |:math:`P(b_n)`                          |
|                                     |                                        |
|[Select old length/angle :math:`b_o`]|:math:`P(b_o)`                          |
+-------------------------------------+----------------------------------------+
|Select new orientation               |:math:`1/\omega`                        |
|                                     |                                        |
|[Select old orientation]             |:math:`1/\omega`                        |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`\min(1, \chi)`                   |
|                                     |                                        |
|[Accept]                             |:math:`[\min(1, 1/\chi)]`               |
+-------------------------------------+----------------------------------------+

Applying (local) detailed balance yields

:math:`\chi = \frac{P(b_o)}{P(b_n)}\exp[-\beta(U_n - U_o)]`.

For example, the probabililty of selecting a bond length/angle, b, by an
entirely random placement in space is given by

:math:`P \propto 1`.

In this case,

:math:`\chi = \exp[-\beta(U_n - U_o)]`

where the energies include all intramolecular and intermolecular terms.
For strong bonds, this becomes highly inefficient because most trials will be
rejected due to the intramolecular term, which is inexpensive to compute,
but the intermolecular term has to be computed for each trial and is expensive.
This is described as scheme one on page 347 of Frenkel and Smit's Understanding
Molecular Simulation.

The second, more efficient scheme is to select the bond/angle by taking the
intramolecular bond energy, :math:`U^b`, into account during the selection.

:math:`P \propto \exp(-\beta U^b)`.

In this case,

:math:`\chi = \frac{P(b_o)}{P(b_n)}\exp\{-\beta[(U_n-U^b_n) - (U_o-U^b_o)]\}`.

Because :math:`U^b` is already taken into account during selection of the
bond/angle using this Rosenbluth form for :math:`P(l)`, :math:`U^b` is
excluded from the energy in the acceptance probability.

\endrst
 */
std::shared_ptr<TrialFactory> MakeTrialGrow(std::vector<argtype> args,
  /// Optionally, set the default values for the following TrialStage arguments:
  /// num_steps, reference_index and new_only.
  /// Any option applied by the above args overwrites this option.
  const argtype& default_args = argtype());

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_GROW_H_

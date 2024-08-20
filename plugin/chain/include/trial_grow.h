#ifndef FEASST_CHAIN_TRIAL_GROW_H_
#define FEASST_CHAIN_TRIAL_GROW_H_

#include <memory>
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

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
  - gibbs_transfer: if true, create two trials which transfer particles between
    configuration_index and configuration_index2 (default: false).
  - regrow: if true, place anywhere in the domain (default: false).
  - transfer_avb: if true, same as transfer but with TrialAddAVB/TrialRemoveAVB
    for the first stage (default: false).
  - regrow_avb2: if true, regrow using TrialAVB2 in the first stage (default: false).
  - regrow_avb4: if true, regrow using TrialAVB4 in the first stage (default: false).
  - translate: if true (default: false), translate site (which is required arg
    for TrialSelectParticle).
    In addition, must have number of stages equal to number of sites.
  - add: if true, create an add trial (default: false).
  - remove: if true, create a remove trial (default: false).
  - add_avb: if true, create an avb add trial (default: false).
  - remove_avb: if true, create an avb remove trial (default: false).
  - default_num_steps: optional default number of steps for all stages.
  - default_reference_index: optional default reference index for all stages.
  - default_new_only: optional default new only for all stages.

  The following options may be used in any argtype.
  If used in the first, then it is a partial regrowth move.
  - bond: if true, add TrialSelectBond and PerturbDistance.
    Requires arguments described in TrialSelectBond.
  - angle: if true, adds TrialSelectAngle and PerturbDistanceAngle.
    Requires arguments described in TrialSelectAngle and TrialSelectBond.
  - dihedral: if true, adds TrialSelectDihedral and PerturbDihedral.
    Requires arguments described in TrialSelectDihedral, TrialSelectAngle
    and TrialSelectBond.
  - branch: if true, adds SelectBranch and PerturbBranch.
    Requires arguments described in SelectBranch, TrialSelectAngle and
    TrialSelectBond.
    Note, only the "2-branch" case, where two sites are placed simultaneous, as
    described in
    https://doi.org/10.1021/acs.jctc.7b00173
    is implemented here.
  - reptate: if true, add TrialSelectBond and PerturbToAnchor.
    Requires arguments described in TrialSelectBond.
  - position_swap: if true, add SelectTwoSites and PerturbPositionSwap.
    Requires arguments described in SelectTwoSites.
  - rigid_body_connector: if true, add TrialSelectBond and PerturbConnector.
    Requires arguments described in TrialSelectBond.
  - rigid_body_angle: if true, add TrialSelectAngle and PerturbDistanceAngleConnector.
    Requires arguments described in TrialSelectAngle.
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
class TrialGrow : public TrialFactoryNamed {
 public:
  TrialGrow() : TrialFactoryNamed() {}
  /// list of arguments, one for each stage.
  TrialGrow(std::vector<argtype> args);
  virtual ~TrialGrow() {}

 protected:
  void build_(std::vector<argtype> * args);
};

inline std::shared_ptr<TrialGrow> MakeTrialGrow(std::vector<argtype> args) {
  return std::make_shared<TrialGrow>(args); }

class TrialGrowFile : public TrialGrow {
 public:
  TrialGrowFile() : TrialGrow() {} // only use for deserialize_map.
  //@{
  /** @name Arguments
    - grow_file: name of TrialGrowFile file with the following format:

      line1: TrialGrowFile

      line2: optional space

      line3: stage with key pair separated by space (e.g., "transfer true site 0")

      lineX: additional stages until end of file or empty line.

      lineY: additional trials with additional stages as described above,
             separated by empty lines.

      See the tutorials for examples.
   */
  explicit TrialGrowFile(argtype args);
  explicit TrialGrowFile(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialGrowFile>(args); }
  virtual ~TrialGrowFile() {}

  //@}
 private:
  void add_(const argtype add_args, std::vector<argtype> * args);
};

inline std::shared_ptr<TrialGrowFile> MakeTrialGrowFile(argtype args) {
  return std::make_shared<TrialGrowFile>(args); }

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_GROW_H_

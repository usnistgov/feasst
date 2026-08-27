
#ifndef FEASST_RXNAVB_TRIAL_RXVB_H_
#define FEASST_RXNAVB_TRIAL_RXVB_H_

#include <memory>
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
\rst

.. image:: ../doc/rxnavb.png

Attempt to react a pair of particles, A and B, with a distance bias.
When A changes to C, B is chosen in the aggregation volume (AV) of A to change
to D, where C and D are given random orientations.
This is expected to improve sampling when particles A and B form attractive
pairs.

A :math:`\rightarrow` C

B(in AV of A) :math:`\rightarrow` D

This Trial supports rigid particles only due to how random orientations are
implemented.

An AV is defined by a spherical shell of given inner and outer radius about a
site in particle A with index of site_index_a.
Particle B is in the AV if site_index_b in particle B is inside the AV
(e.g., these are site-based definitions of AV).

Because the particle types change, this Trial is derived in the semi-grand
canonical ensemble.
The derivation of the acceptance criteria is dervied by a similar method as
described in https://github.com/usnistgov/best-practices-mc:

+-------------------------------------+----------------------------------------+
|Forward event                        |:math:`\alpha_{o\rightarrow n}`         |
+-------------------------------------+----------------------------------------+
|Choose A and B reaction              |:math:`1/2`                             |
+-------------------------------------+----------------------------------------+
|Choose particle of type A            |:math:`1/N_A`                           |
+-------------------------------------+----------------------------------------+
|Choose site_index_b in AV of         |:math:`1/N^{AV}_B`                      |
|chosen site_index_a                  |                                        |
+-------------------------------------+----------------------------------------+
|Change particle A to C and reorient  |:math:`P_{\omega C}`                    |
+-------------------------------------+----------------------------------------+
|Change particle B to D and reorient  |:math:`P_{\omega D}`                    |
+-------------------------------------+----------------------------------------+

+-------------------------------------+----------------------------------------+
|Reverse event                        |:math:`\alpha_{n \rightarrow o}`        |
+-------------------------------------+----------------------------------------+
|Choose C and D reaction              |:math:`1/2`                             |
+-------------------------------------+----------------------------------------+
|Choose particle of type C            |:math:`1/(N_C+1)`                       |
+-------------------------------------+----------------------------------------+
|Choose site_index_d in AV of         |:math:`1/(N^{AV}_D+1)`                  |
|chosen site_index_c                  |                                        |
+-------------------------------------+----------------------------------------+
|Change particle C to A and reorient  |:math:`P_{\omega A}`                    |
+-------------------------------------+----------------------------------------+
|Change particle D to B and reorient  |:math:`P_{\omega B}`                    |
+-------------------------------------+----------------------------------------+

Application of local detailed balance yields the acceptance probability,

:math:`\chi = \frac{P_{\omega A} P_{\omega B}}{P_{\omega C} P_{\omega D}}
\frac{z_C z_D}{z_A z_B} \frac{N_A N^{AV}_B}{(N_C + 1)(N^{AV}_D+1)}
e^{-\beta\Delta U}`

where :math:`N^{AV}_D` is the number of site_index_d of particle type D in the
AV in the old configuration (before reaction).

The reverse reaction reaction is implemented as follows:

C :math:`\rightarrow` A

D(in AV of C) :math:`\rightarrow` B

The acceptance of this reverse move is

:math:`\chi = \frac{P_{\omega C} P_{\omega D}}{P_{\omega A} P_{\omega B}}
\frac{z_A z_B}{z_C z_D} \frac{N_C N^{AV}_D}{(N_A + 1)(N^{AV}_B + 1)}
e^{-\beta\Delta U}`
\endrst

In order to obey detailed balance, the AV of A and C must be identical.
Otherwise, the reverse trial may not be possible.

The current implementation assumes the AVB site is the first on each particle.
*/
class TrialRxVB : public TrialFactoryNamed {
 public:
  //@{
  /** @name Arguments
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
    - Trial arguments.
    - TrialStage arguments.
    - SelectParticleAVB arguments
      - target_particle_type and target_site defines the AV
      - particle_type and site (from TrialSelectParticle) defines particle in AV.
    - target_particle_type_morph: type of particle to change target
    - particle_type_morph: type of particle to change particle in AV
   */
  explicit TrialRxVB(argtype args = argtype());
  explicit TrialRxVB(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialRxVB>(args); }
  virtual ~TrialRxVB() {}
  //@}
};

/**
  Attempt either A + B(in AV of A) -> C + D
  or C + D(in AV of C) -> A + B
 */
class TrialRxVBHalf : public Trial {
 public:
  /**
    args:
    - avb: if true, use out->in. Otherwise, in->out.
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
   */
  explicit TrialRxVBHalf(argtype args = argtype());
  explicit TrialRxVBHalf(argtype * args);

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialRxVBHalf>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialRxVBHalf>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialRxVBHalf(std::istream& istr);
  virtual ~TrialRxVBHalf() {}
};

// Process args, which can also be used in TrialGrow
void gen_rxnavb_args_(const bool avb, argtype * args, argtype * perturb_args);

}  // namespace feasst

#endif  // FEASST_RXNAVB_TRIAL_RXVB_H_

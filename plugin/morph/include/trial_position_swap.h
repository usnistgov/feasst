
#ifndef FEASST_MORPH_TRIAL_POSITION_SWAP_H_
#define FEASST_MORPH_TRIAL_POSITION_SWAP_H_

#include <memory>
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
Attempt to change the identity of two particles of different types in
the same configuration.
This trial is only valid for rigid particles.
This implementation is extremely similar to TrialPositionSwap over different
configurations.

\rst
The limiting distribution in the canonical ensemble is given by

:math:`\frac{\Pi_{n}}{\Pi_{o}} = e^{-\beta\Delta (U_1+U_2)}`

See (https://doi.org/10.33011/livecoms.6.1.3289)

The transition probabilities are as follows.

+-------------------------------------+---------------------------------------------------+
|Forward                              |:math:`\pi_{o \rightarrow n}`                      |
+-------------------------------------+---------------------------------------------------+
|Choose particle of type i            |:math:`1/N_{i}`                                    |
+-------------------------------------+---------------------------------------------------+
|Choose particle of type j            |:math:`1/N_{j}`                                    |
+-------------------------------------+---------------------------------------------------+
|Exchange the types of the particles  |:math:`P_{\omega i n}\mathrm{d}\boldsymbol{\omega}`|
|and rotate both of them              |:math:`P_{\omega j n}\mathrm{d}\boldsymbol{\omega}`|
+-------------------------------------+---------------------------------------------------+

+-------------------------------------+---------------------------------------------------+
|Reverse                              |:math:`\pi_{n \rightarrow o}`                      |
+-------------------------------------+---------------------------------------------------+
|Choose particle of type i            |:math:`1/N_{i}`                                    |
+-------------------------------------+---------------------------------------------------+
|Choose particle of type j            |:math:`1/N_{j}`                                    |
+-------------------------------------+---------------------------------------------------+
|Exchange the types of the particles  |:math:`P_{\omega i o}\mathrm{d}\boldsymbol{\omega}`|
|and rotate both of them              |:math:`P_{\omega j o}\mathrm{d}\boldsymbol{\omega}`|
+-------------------------------------+---------------------------------------------------+

Application of local detailed balance yields the acceptance probability,

:math:`\chi = \frac{P_{\omega i o}P_{\omega j o}}{P_{\omega i n}P_{\omega j n}} e^{-\beta \Delta U}`
\endrst
 */
class TrialPositionSwap : public Trial {
 public:
  /**
    args:
    - TrialSelectParticle arguments.
    - particle_type_morph: the name of the particle type to morph into (and
      vice versa), not equal to particle_type.
   */
  explicit TrialPositionSwap(argtype args = argtype());
  explicit TrialPositionSwap(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialPositionSwap>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialPositionSwap>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialPositionSwap(std::istream& istr);
  virtual ~TrialPositionSwap() {}
};

}  // namespace feasst

#endif  // FEASST_MORPH_TRIAL_POSITION_SWAP_H_

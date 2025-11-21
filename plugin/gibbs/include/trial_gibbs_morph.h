
#ifndef FEASST_GIBBS_TRIAL_GIBBS_MORPH_H_
#define FEASST_GIBBS_TRIAL_GIBBS_MORPH_H_

#include <memory>
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
Attempt to change the identity of two particles of different types in
different configurations.
This trial is only valid for rigid particles.

\rst
The limiting distribution in the Gibbs ensemble with changing number of
particles is given by

:math:`\frac{\Pi_{n}}{\Pi_{o}} = V_1^{N_{ni1}-N_{oi1}} V_2^{N_{ni2}-N_{oi2}}e^{-\beta\Delta (U_1+U_2)}`

See (https://doi.org/10.33011/livecoms.6.1.3289)

The transition probabilities are as follows.

+-------------------------------------+--------------------------------------------------+
|Forward                              |:math:`\pi_{o \rightarrow n}`                     |
+-------------------------------------+--------------------------------------------------+
|Choose particle of type i in domain 1|:math:`1/N_{i1}`                                  |
+-------------------------------------+--------------------------------------------------+
|Choose particle of type j in domain 2|:math:`1/N_{j2}`                                  |
+-------------------------------------+--------------------------------------------------+
|Exchange the types of the particles  |:math:`P_{\omega i2}\mathrm{d}\boldsymbol{\omega}`|
|and rotate both of them              |:math:`P_{\omega j1}\mathrm{d}\boldsymbol{\omega}`|
+-------------------------------------+--------------------------------------------------+

+-------------------------------------+--------------------------------------------------+
|Reverse                              |:math:`\pi_{n \rightarrow o}`                     |
+-------------------------------------+--------------------------------------------------+
|Choose particle of type i in domain 2|:math:`1/(N_{i2}+1)`                              |
+-------------------------------------+--------------------------------------------------+
|Choose particle of type j in domain 1|:math:`1/(N_{j1}+1)`                              |
+-------------------------------------+--------------------------------------------------+
|Exchange the types of the particles  |:math:`P_{\omega i1}\mathrm{d}\boldsymbol{\omega}`|
|and rotate both of them              |:math:`P_{\omega j2}\mathrm{d}\boldsymbol{\omega}`|
+-------------------------------------+--------------------------------------------------+

Application of local detailed balance yields the acceptance probability,

:math:`\chi = \frac{N_{i1} N_{j2}}{(N_{i2}+1)(N_{j1}+1)} \frac{P_{\omega i1}P_{\omega j2}}{P_{\omega i2}P_{\omega j1}} e^{-\beta (\Delta U_1 + \Delta U_2)}`
\endrst

Attempt TrialGibbsMorphOneWay with equal probability in either direction.
 */
class TrialGibbsMorph : public TrialFactoryNamed {
 public:
  //@{
  /** @name Arguments
    - configs: two comma-separated names of configurations to transfer between
      (default: "0,1").
    - TrialSelectParticle arguments.
    - particle_type_morph: the name of the particle type to morph into (and
      vice versa), not equal to particle_type.
   */
  explicit TrialGibbsMorph(argtype args = argtype());
  explicit TrialGibbsMorph(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialGibbsMorph>(args); }
  virtual ~TrialGibbsMorph() {}
  //@}
};

/// Attempt to transfer a particle from one configuration to another.
class TrialGibbsMorphOneWay : public Trial {
 public:
  /**
    args:
    - to_config: name of configuration to send the particle.
    - config: from TrialSelect, configuration which donates a
      particle (default: 0).
   */
  explicit TrialGibbsMorphOneWay(argtype args = argtype());
  explicit TrialGibbsMorphOneWay(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialGibbsMorphOneWay>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialGibbsMorphOneWay>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialGibbsMorphOneWay(std::istream& istr);
  virtual ~TrialGibbsMorphOneWay() {}
};

}  // namespace feasst

#endif  // FEASST_GIBBS_TRIAL_GIBBS_MORPH_H_

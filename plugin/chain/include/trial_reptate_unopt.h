#ifndef FEASST_CHAIN_TRIAL_REPTATE_UNOPTH_
#define FEASST_CHAIN_TRIAL_REPTATE_UNOPTH_

#include <memory>
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
\rst

Attempt to reptate a linear chain particle of a given type as follows.
Randomly pick one of two ends of the chain.
Starting from that end, move each site to the next, and regrow the position of
the last.
This reptation trial is general to most linear particles and does not use
optimization which assumes that all beads on the chain are identical
(See TrialReptate).

+-------------------------------------+----------------------------------------+
|Forward event                        |:math:`\alpha_{o\rightarrow n}`         |
+-------------------------------------+----------------------------------------+
|Choose particle of type A            |:math:`1/N_A`                           |
+-------------------------------------+----------------------------------------+
|Choose end of linear particle        |:math:`1/2`                             |
+-------------------------------------+----------------------------------------+
|Move each site to the next and pick a|:math:`P_{\omega n}`                    |
|new orientation for the last site.   |                                        |
+-------------------------------------+----------------------------------------+

+-------------------------------------+----------------------------------------+
|Reverse event                        |:math:`\alpha_{n\rightarrow o}`         |
+-------------------------------------+----------------------------------------+
|Choose particle of type A            |:math:`1/N_A`                           |
+-------------------------------------+----------------------------------------+
|Choose opposite end of linear        |:math:`1/2`                             |
|particle                             |                                        |
+-------------------------------------+----------------------------------------+
|Move each site to the next and pick a|:math:`P_{\omega o}`                    |
|new orientation for the site on the  |                                        |
|opposite end                         |                                        |
+-------------------------------------+----------------------------------------+

Because this trial only involves changing particle positions, the canonical
ensemble microstate probability distribution is utilized in the application of
local detailed balance.

:math:`\chi = \frac{P_{\omega o}}{P_{\omega n}}e^{-\beta\Delta U}`

\endrst
*/
class TrialReptateUnopt : public TrialFactoryNamed {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - TrialStage arguments.
    - particle_type: See TrialSelect.
    - config: See TrialSelect.
   */
  explicit TrialReptateUnopt(argtype args = argtype());
  explicit TrialReptateUnopt(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialReptateUnopt>(args); }
  virtual ~TrialReptateUnopt() {}
  //@}
};

class TrialReptateUnoptHalf : public Trial {
 public:
  /**
    args:
    - reverse: if 1, reverse direction (default: 0).
   */
  explicit TrialReptateUnoptHalf(argtype args = argtype());
  explicit TrialReptateUnoptHalf(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute(Criteria * criteria, System * system) override;

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialReptateUnoptHalf>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialReptateUnoptHalf>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialReptateUnoptHalf(std::istream& istr);
  virtual ~TrialReptateUnoptHalf();
  //@}
 private:
  argtype args_;
};

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_REPTATE_UNOPTH_

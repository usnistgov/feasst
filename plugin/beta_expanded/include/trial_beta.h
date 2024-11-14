#ifndef FEASST_BETA_EXPANDED_TRIAL_BETA_H_
#define FEASST_BETA_EXPANDED_TRIAL_BETA_H_

#include <memory>
#include "monte_carlo/include/trial.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
Attempt to change the inverse temperature, \f$\beta=\frac{1}{k_B T}\f$ by a
fixed amount.
\rst

The derivation of the acceptance criteria follows a similar procedure as
descibed in TrialComputeMove, except with the following differences.

The limiting distribution in the canonical canonical ensemble is

:math:`\pi_i \propto e^{-\beta U_i}`

The transition probabilities are as follows.

+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Select increase of decrease          |:math:`1/2`                             |
|                                     |                                        |
|[Select increase of decrease]        |:math:`[1/2]`                           |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

Application of local detailed balance yields the acceptance probability,
:math:`\chi`.

:math:`\frac{e^{-\beta_o U}}{2}min(1, \chi) =
\frac{e^{-\beta_n U}}{2} min(1, 1/\chi)`

:math:`\chi = e^{-(\beta_n - \beta_o) U}`

\endrst
 */
class TrialBeta : public Trial {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - PerturbBeta arguments.
    - Tunable arguments.
   */
  explicit TrialBeta(argtype args = argtype());
  explicit TrialBeta(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialBeta>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialBeta>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialBeta(std::istream& istr);
  virtual ~TrialBeta() {}
  //@}
};

}  // namespace feasst

#endif  // FEASST_BETA_EXPANDED_TRIAL_BETA_H_


#ifndef FEASST_BETA_EXPANDED_COMPUTE_BETA_H_
#define FEASST_BETA_EXPANDED_COMPUTE_BETA_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

/**
\rst
Attempt to change the inverse temperature, :math:`\beta`.

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

:math:`\chi = \frac{e^{-(\beta_n - \beta_o) U}}`

\endrst
 */
class ComputeBeta : public TrialCompute {
 public:
  ComputeBeta();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeBeta(std::istream& istr);
  virtual ~ComputeBeta() {}

 protected:
  void serialize_trial_compute_add_(std::ostream& ostr) const;
};

inline std::shared_ptr<ComputeBeta> MakeComputeBeta() {
  return std::make_shared<ComputeBeta>();
}
}  // namespace feasst

#endif  // FEASST_BETA_EXPANDED_COMPUTE_BETA_H_

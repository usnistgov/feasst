
#ifndef FEASST_MODEL_EXPANDED_COMPUTE_MODEL_H_
#define FEASST_MODEL_EXPANDED_COMPUTE_MODEL_H_

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
Attempt to change the Model::model_index.

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

:math:`\frac{e^{-\beta U_0}}{2}min(1, \chi) =
\frac{e^{-\beta U_n}}{2} min(1, 1/\chi)`

:math:`\chi = e^{-\beta(U_n - U_o)}`

\endrst
 */
class ComputeModel : public TrialCompute {
 public:
  ComputeModel();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeModel(std::istream& istr);
  virtual ~ComputeModel() {}

 protected:
  void serialize_trial_compute_add_(std::ostream& ostr) const;
};

inline std::shared_ptr<ComputeModel> MakeComputeModel() {
  return std::make_shared<ComputeModel>();
}
}  // namespace feasst

#endif  // FEASST_MODEL_EXPANDED_COMPUTE_MODEL_H_

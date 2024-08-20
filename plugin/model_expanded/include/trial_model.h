#ifndef FEASST_MODEL_EXPANDED_TRIAL_MODEL_H_
#define FEASST_MODEL_EXPANDED_TRIAL_MODEL_H_

#include <memory>
#include "monte_carlo/include/trial.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
\rst
Attempt to change the Model::model_index by plus or minus one.

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
class TrialModel : public Trial {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - PerturbModel arguments.
    - Tunable arguments.
   */
  explicit TrialModel(argtype args = argtype());
  explicit TrialModel(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialModel>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialModel>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialModel(std::istream& istr);
  virtual ~TrialModel() {}
  //@}
};

inline std::shared_ptr<TrialModel> MakeTrialModel(argtype args = argtype()) {
  return std::make_shared<TrialModel>(args); }

}  // namespace feasst

#endif  // FEASST_MODEL_EXPANDED_TRIAL_MODEL_H_

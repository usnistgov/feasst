
#ifndef FEASST_MONTE_CARLO_COMPUTE_ADD_MULTIPLE_H_
#define FEASST_MONTE_CARLO_COMPUTE_ADD_MULTIPLE_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

/**
Attempt to add multiple particles.

For a derivation of the acceptance criteria, see TrialComputeMove and
TrialComputeAdd for reference.
For m particles added of types t1, t2, ..., tm:

\rst
+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Select insert trial                  |:math:`1/w`                             |
|                                     |                                        |
|[Select remove trial]                |:math:`[1/w]`                           |
+-------------------------------------+----------------------------------------+
|Place particle of type t1            |:math:`1/V`                             |
|                                     |                                        |
|[Delete particle type t1]            |:math:`\left[\frac{1}{N_{t1}+1}\right]` |
+-------------------------------------+----------------------------------------+
|Place particle of type t2            |:math:`1/V`                             |
|                                     |                                        |
|[Delete particle type t2]            |:math:`\left[\frac{1}{N_{t2}+1}\right]` |
+-------------------------------------+----------------------------------------+
| ...                                 | ...                                    |
+-------------------------------------+----------------------------------------+
|Place particle of type tm            |:math:`1/V`                             |
|                                     |                                        |
|[Delete particle type tm]            |:math:`\left[\frac{1}{N_{tm}+1}\right]` |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

Application of local detailed balance yeilds the acceptance probability,

:math:`\chi = e^{-\beta\Delta U}\prod_{i=1,...,m}\frac{Ve^{\beta\mu_i}}
{(N_{ti}+1)\Lambda^d}`
\endrst
 */
class ComputeAddMultiple : public TrialCompute {
 public:
  ComputeAddMultiple();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeAddMultiple(std::istream& istr);
  virtual ~ComputeAddMultiple() {}

 protected:
  void serialize_compute_add_multiple_(std::ostream& ostr) const;
};

inline std::shared_ptr<ComputeAddMultiple> MakeComputeAddMultiple() {
  return std::make_shared<ComputeAddMultiple>();
}
}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_COMPUTE_ADD_MULTIPLE_H_

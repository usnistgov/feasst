
#ifndef FEASST_CHARGE_COMPUTE_ADD_MULTIPLE_H_
#define FEASST_CHARGE_COMPUTE_ADD_MULTIPLE_H_

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
For adding multiple particles i, j, ..., z:

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
|Place particle of type i             |:math:`1/V`                             |
|                                     |                                        |
|[Delete particle type i]             |:math:`\left[1/(N_{i}+                  |
|                                     |\sum_{a=i}^{z}\delta_{ia})\right]`      |
+-------------------------------------+----------------------------------------+
|Place particle of type j             |:math:`1/V`                             |
|                                     |                                        |
|[Delete particle type j]             |:math:`\left[1/(N_{j}+                  |
|                                     |\sum_{a=j}^{z}\delta_{ja})\right]`      |
+-------------------------------------+----------------------------------------+
| ...                                 | ...                                    |
+-------------------------------------+----------------------------------------+
|Place particle of type z             |:math:`1/V`                             |
|                                     |                                        |
|[Delete particle type z]             |:math:`\left[1/(N_{z}+1)\right]`        |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

where :math:`\delta_{ab} = 1` when the types of particles :math:`a` and
:math:`b` are identical.
Otherwise, :math:`\delta = 0`.

Application of local detailed balance yields the acceptance probability,

:math:`\chi = e^{-\beta\Delta U}\prod_{a=i}^z\frac{Ve^{\beta\mu_a}}
{(N_a+\sum_{b=a}^{z}\delta_{ab})\Lambda^d}`

This equation was derived from the perspective of the old state.
If the new state is the perspective,
:math:`N_a+\sum_{b=a}^{z}\delta_{ab} \rightarrow N_a-\sum_{b=i}^{a-1}\delta_{ab}`

\endrst
 */
class ComputeAddMultiple : public TrialCompute {
 public:
  /**
    args:
    - shift: macrostate shift (default: -1).
   */
  explicit ComputeAddMultiple(argtype args = argtype());
  explicit ComputeAddMultiple(argtype * args);

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
  int shift_;

 private:
  std::vector<int> delta_;
};

inline std::shared_ptr<ComputeAddMultiple> MakeComputeAddMultiple() {
  return std::make_shared<ComputeAddMultiple>();
}
}  // namespace feasst

#endif  // FEASST_CHARGE_COMPUTE_ADD_MULTIPLE_H_

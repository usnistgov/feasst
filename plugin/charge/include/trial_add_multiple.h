
#ifndef FEASST_CHARGE_TRIAL_ADD_MULTIPLE_H_
#define FEASST_CHARGE_TRIAL_ADD_MULTIPLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {


// parse the number of particle types.
std::vector<int> ptypes(argtype * args);

/**
Attempt to add multiple particles.
Typically requires the use of a reference index.

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
class TrialAddMultiple : public Trial {
 public:
  //@{
  /** @name Arguments
    - particle_type[i]: the i-th type of particle to add.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
    - TrialStage arguments.
   */
  explicit TrialAddMultiple(argtype args = argtype());
  explicit TrialAddMultiple(argtype * args);
  //@}
  /** @name Public Functions
   */
  //@{
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAddMultiple>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAddMultiple>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAddMultiple(std::istream& istr);
  virtual ~TrialAddMultiple() {}
  //@}
};

inline std::shared_ptr<TrialAddMultiple> MakeTrialAddMultiple(argtype args = argtype()) {
  return std::make_shared<TrialAddMultiple>(args); }

}  // namespace feasst

#endif  // FEASST_CHARGE_TRIAL_ADD_MULTIPLE_H_

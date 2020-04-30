
#ifndef FEASST_MONTE_CARLO_TRIAL_COMPUTE_MOVE_H_
#define FEASST_MONTE_CARLO_TRIAL_COMPUTE_MOVE_H_

#include <vector>
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

/**
Move a selection of particles and sites.

Derivation of the acceptance criteria follows Lecture 9 of Prof. David
Kofke's Molecular Simulation course CE 530.

http://www.eng.buffalo.edu/~kofke/ce530/Lectures/lectures.html

For this type of move, the potential energy of the system, U, is the only
thermodynamic variable which changes.
Thus, whether in the canonical ensemble or otherwise, the probability
distribution,

\rst
:math:`\pi_i \propto e^{-\beta U_i}`

because the other thermodynamic variables such as number of particles,
or volume do not change.

The following table describes the transition probabilities associated with
the chosen trial move, and its reverse trial that is considered for the
purpose of satisfying detailed balance.
Thus, the probability shown represents the probability of transitioning from
the old to the new state, :math:`\pi_{on}`.
And the reverse transition probability is from new to old, :math:`\pi_{no}`.

+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Select molecule of type t            |:math:`1/N_t`                           |
|                                     |                                        |
|[Select molecule of type t]          |:math:`[1/N_t]`                         |
+-------------------------------------+----------------------------------------+
|Move to new position                 |:math:`1/v`                             |
|                                     |                                        |
|[Move to old position]               |:math:`[1/v]`                           |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

where :math:`\chi` is the acceptance probability that we can now derive by
imposing the (local) detailed balance condition.

For (local) detailed balance, the probability of being in the old state,
:math:`\pi_o`, times the probability of transitioning from the old to the new
state, :math:`\pi_{on}`, should be equal to the probability of being in the new
state, :math:`\pi_n` times the probability of transitioning from the new to old
state, :math:`\pi_{no}`.

:math:`\pi_o \pi_{on} = \pi_n \pi_{no}`

Substituting the probability distribution and transition probabilities yeilds

:math:`\frac{e^{-\beta U_o}min(1, \chi)}{N_t v} =
\frac{e^{-\beta U_n}min(1, 1/\chi)}{N_t v}`

:math:`min(1, \chi)/min(1, 1/\chi) = e^{-\beta(U_n - U_o)} =
e^{-\beta\Delta U}`

The left hand side is :math:`\chi` for both cases of :math:`\chi <= 1`
and :math:`\chi > 1`. Thus,

:math:`\chi = e^{-\beta\Delta U}`
\endrst
 */
class TrialComputeMove : public TrialCompute {
 public:
  TrialComputeMove();

  void perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) override;
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeMove(std::istream& istr);
  virtual ~TrialComputeMove() {}

 protected:
  void serialize_trial_compute_move_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialComputeMove> MakeTrialComputeMove() {
  return std::make_shared<TrialComputeMove>();
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_COMPUTE_MOVE_H_

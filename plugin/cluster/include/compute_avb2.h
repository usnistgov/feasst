
#ifndef FEASST_MONTE_CARLO_COMPUTE_AVB2_H_
#define FEASST_MONTE_CARLO_COMPUTE_AVB2_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute_move.h"

namespace feasst {

/**
This trial move is refered to as AVB2 as described in the following manuscript:

https://doi.org/10.1021/jp012209k

First, a target particle of type "t" is selected with a site with index
site_index_t.
With probability, \f$P_{bias}\f$, a particle of type "a" with site_index_a
outside the aggregation volume (AV) of site_index_t is placed inside the AV.
Otherwise, site_index_a inside the AV is placed outisde the AV.

The derivation of the acceptance probability follows that as described in both
TrialComputeMove and ComputeAddAVB with the following transition probabilities.

For this type of move, the potential energy of the system, U, is the only
thermodynamic variable which changes.
Thus, whether in the canonical ensemble or otherwise, the probability
distribution,

\rst
:math:`\pi_i \propto e^{-\beta U_i}`

+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Select particle of type t            |:math:`1/N_t`                           |
|                                     |                                        |
|[Select particle of type t]          |:math:`[1/N_t]`                         |
+-------------------------------------+----------------------------------------+
|Select out->in move                  |:math:`P_{bias}`                        |
|                                     |                                        |
|[Select in->out move]                |:math:`1 - P_{bias}`                    |
+-------------------------------------+----------------------------------------+
|Select site_index_a not in AV of     |:math:`1/(N_a - N^{s,AV}_a -            |
|site_index_t                         |\delta_{ta})`                           |
|                                     |                                        |
|[Select site_index_a in AV of        |:math:`1/(N^{s,AV}_a + 1)`              |
|site_index_t]                        |                                        |
+-------------------------------------+----------------------------------------+
|Place site_index_a inside AV of      |:math:`1/v_{AV}`                        |
|site_index_t                         |                                        |
|                                     |                                        |
|[Place site_index_a outside AV of    |:math:`1/(V - v_{AV})`                  |
|site_index_t]                        |                                        |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

where :math:`N_t` is the number of particles of type t,
:math:`v_{AV}` is the aggregation volume,
:math:`N^{s,AV}_a` is the number of sites with site_index_a, in particles of
type a, that are in the AV of site_index_t,
:math:`V` is the volume of the entire system
and :math:`\delta_{ta}` is the Kronecker delta as a function of the type index
of the target and added particles (i.e., if add and target are of same type,
then :math:`\delta_{ta} = 1`, otherwise :math:`\delta_{ta} = 0`).

Application of local detailed balance yields the following acceptance
probability for the out->in move:

:math:`\chi =
\frac{1-P_{bias}}{P_{bias}}
\frac{N_a - N^{s,AV}_a - \delta_{ta}}{N^{s,AV}_a + 1}
\frac{v_{AV}}{V - v_{AV}}
e^{-\beta\Delta U}`.

The in->out move can be derived by the switching the foward and reverse moves.
In that case, one is no longer added to :math:`N^{s,AV}_a`, but to
:math:`N_a - N^{s,AV}_a - \delta_{ta}` instead.
The resulting acceptance probability is

:math:`\chi =
\frac{P_{bias}}{1-P_{bias}}
\frac{N^{s,AV}_a}{N_a - N^{s,AV}_a - \delta_{ta} + 1}
\frac{V - v_{AV}}{v_{AV}}
e^{-\beta\Delta U}`.

\endrst
 */
class ComputeAVB2 : public TrialComputeMove {
 public:
  /**
    args:
    - probability_out_to_in: probability to attempt an out->in move
      (default: 0.5).
    - out_to_in: true if out->in move, otherwise false.
   */
  explicit ComputeAVB2(argtype args = argtype());
  explicit ComputeAVB2(argtype * args);

  /// Return the probability of an out to in move.
  double probability_out_to_in() const { return p_bias_; }

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeAVB2(std::istream& istr);
  virtual ~ComputeAVB2() {}

 protected:
  void serialize_compute_avb2_(std::ostream& ostr) const;
  double p_bias_;
  bool out_to_in_;
};

inline std::shared_ptr<ComputeAVB2> MakeComputeAVB2(argtype args = argtype()) {
  return std::make_shared<ComputeAVB2>(args);
}
}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_COMPUTE_AVB_H_

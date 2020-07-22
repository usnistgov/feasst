
#ifndef FEASST_CLUSTER_COMPUTE_ADD_AVB_DIVALENT_H_
#define FEASST_CLUSTER_COMPUTE_ADD_AVB_DIVALENT_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {


/**
Attempt to add a particle of type "t" with site_index_t anywhere in the Domain.
Then add a second particle of type "a" with site_index_a in the AV of
site_index_t.
Finally, add a third particle of type "b" with site_index_a in the AV of
site_index_t.

The derivation of the acceptance criteria follows a similar procedure as
described in ComputeAddAVB.

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
|Insert site_index_t in box           |:math:`1/V`                             |
|                                     |                                        |
|[Remove site_index_t]                |:math:`1/(N_t + 1 + \delta_{ta} +       |
|                                     |\delta_{tb})`                           |
+-------------------------------------+----------------------------------------+
|Insert site_index_a in AV of         |:math:`1/v_{AV}`                        |
|site_index_t                         |                                        |
|                                     |                                        |
|[Remove site_index_a in AV of        |:math:`1/(N^{s,AV}_a + 1 + \delta_{ab})`|
|site_index_t                         |                                        |
+-------------------------------------+----------------------------------------+
|Insert site_index_b in AV of         |:math:`1/v_{AV}`                        |
|site_index_t                         |                                        |
|                                     |                                        |
|[Remove site_index_b in the AV of    |:math:`1/(N^{s,AV}_b + 1)`              |
|site_index_t]                        |                                        |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

where :math:`N_t` is the number of particles of type t,
:math:`v_{AV}` is the aggregation volume,
:math:`N^{s,AV}_a` is the number of sites with site_index_a, in particles of
type a, that are in the AV of site_index_t and
:math:`\delta_{ta}` is the Kronecker delta as a function of the type index
of the target and added particles (i.e., if add and target are of same type,
then :math:`\delta_{ta} = 1`, otherwise :math:`\delta_{ta} = 0`).

Application of local detailed balance yields the acceptance probability

:math:`\frac{e^{-\beta U_o + \beta\mu_a N_a + \beta\mu_t N_t}}{\Lambda^{dN} w V v_{AV}^2}
min(1, \chi) =
\frac{e^{-\beta U_n + \beta\mu_a(N_a+2) + \beta\mu_t (N_t+1)}}{\Lambda^{d(N+3)}w(N_t + 1 + \delta_{ta} + \delta_{tb})
(N^{s,AV}_a+1+\delta_{ab})(N^{s,AV}_a+1)}min(1, 1/\chi)`

:math:`\chi =
\frac{V}{N_t + 1 + \delta_{ta} + \delta_{tb}}
\frac{v_{AV}}{N^{s,AV}_a + 1 + \delta_{ab}}
\frac{v_{AV}}{N^{s,AV}_a + 1}
\frac{1}{\Lambda^3d}
e^{-\beta\Delta U + 2\beta\mu_a + \beta\mu_t}`

\endrst
 */
class ComputeAddAVBDivalent : public TrialCompute {
 public:
  explicit ComputeAddAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria);

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeAddAVBDivalent(std::istream& istr);
  virtual ~ComputeAddAVBDivalent() {}

 protected:
  void serialize_compute_add_avb_divalent_(std::ostream& ostr) const;

 private:
  std::shared_ptr<NeighborCriteria> neighbor_criteria_;

  // temporary
  Select neighbors_;
};

inline std::shared_ptr<ComputeAddAVBDivalent> MakeComputeAddAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria) {
  return std::make_shared<ComputeAddAVBDivalent>(neighbor_criteria);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_COMPUTE_ADD_AVB_DIVALENT_H_

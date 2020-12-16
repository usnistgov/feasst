
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

:math:`\chi =
\frac{V}{N_t + 1 + \delta_{ta} + \delta_{tb}}
\frac{v_{AV}}{N^{s,AV}_a + 1 + \delta_{ab}}
\frac{v_{AV}}{N^{s,AV}_a + 1}
\frac{1}{\Lambda^{3d}}
e^{-\beta\Delta U + \beta\mu_a + \beta\mu_b + \beta\mu_t}`

\endrst
 */
class ComputeAddAVBDivalent : public TrialCompute {
 public:
  /**
    args:
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
   */
  explicit ComputeAddAVBDivalent(const argtype& args = argtype());

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
  int neighbor_;

  // temporary
  Select neighbors_;
};

inline std::shared_ptr<ComputeAddAVBDivalent> MakeComputeAddAVBDivalent(
    const argtype& args = argtype()) {
  return std::make_shared<ComputeAddAVBDivalent>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_COMPUTE_ADD_AVB_DIVALENT_H_

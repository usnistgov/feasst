
#ifndef FEASST_CLUSTER_TRIAL_ADD_AVB_H_
#define FEASST_CLUSTER_TRIAL_ADD_AVB_H_

#include <memory>
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
Attempt to add a particle of type "a" (for add) with a bias for site_index_a
to be in the aggregation volume (AV) of site_index_t ("t" for target) of a
particle of type "t".

The derivation of the acceptance probability follows that as described in both
TrialComputeMove and TrialComputeAdd with the following transition
probabilities.

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
|Select particle of type t            |:math:`1/N_t`                           |
|                                     |                                        |
|[Select particle of type t]          |:math:`\left[\frac{1}{N_t+\delta_{ta}}  |
|                                     |\right]`                                |
+-------------------------------------+----------------------------------------+
|Place site_index_a in AV of          |:math:`1/v_{AV}`                        |
|site_index_t                         |                                        |
|                                     |                                        |
|[Remove site_index_a in AV of        |:math:`[1/N^{s,AV}_a + 1]`              |
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

:math:`\frac{e^{-\beta U_o + \beta\mu_a N_a}}{\Lambda^{dN} w N_t v_{AV}}
min(1, \chi) =
\frac{e^{-\beta U_n + \beta\mu_a(N_a+1)}}{\Lambda^{d(N+1)}w(N_t + \delta_{ta})
(N^{s,AV}_a+1)}min(1, 1/\chi)`

:math:`\chi =
\frac{N_t}{N_t+\delta_{ta}}
\frac{v_{AV}}{(N^{s,AV}_a + 1)\Lambda^d}
e^{-\beta\Delta U + \beta\mu_a}`

Note that if the add was already performed, and :math:`\delta_{ta}=1` then
the equations are altered from the perspective of the old to the new state.

:math:`\frac{N_t}{N_t + \delta_{ta}} \rightarrow \frac{N_t - \delta_{ta}}{N_t}`

\endrst
 */
class TrialAddAVB : public Trial {
 public:
  explicit TrialAddAVB(argtype args = argtype());
  explicit TrialAddAVB(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAddAVB>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAddAVB>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAddAVB(std::istream& istr);
  virtual ~TrialAddAVB() {}
};

inline std::shared_ptr<TrialAddAVB> MakeTrialAddAVB(argtype args = argtype()) {
  return std::make_shared<TrialAddAVB>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_ADD_AVB_H_

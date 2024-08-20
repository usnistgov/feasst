
#ifndef FEASST_CLUSTER_TRIAL_AVB4_H_
#define FEASST_CLUSTER_TRIAL_AVB4_H_

#include <memory>
#include "monte_carlo/include/trial.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
AVB4 is an extension of AVB3 as described in the following manuscript:

https://doi.org/10.1021/jp012209k

\rst
But in AVB4, in :math:`\rightarrow` out and out
:math:`\rightarrow` in moves are not considered.
This is because AVB2 handles these moves more efficiently without having to
select a third particle.
One may also allow the :math:`j,k` particles to be either the same or
to have overlapping aggregation volumes (AV).

The AVB4 algorithm proceeds as follow.
First, a particle of type "j" is selected with a site with index site_index_j.
Second, a particle of type "k" is selected with a site with index site_index_k.
With probability, :math:`P_{bias}`, select site_index_a in the AV of
site_index_k and place it in the AV of site_index_j.
Otherwise, select site_index_a in the AV of site_index_j and place it in the AV
of site_index_k.

For this type of move, the potential energy of the system, U, is the only
thermodynamic variable which changes.
Thus, whether in the canonical ensemble or otherwise, the probability
distribution,

:math:`\pi_i \propto e^{-\beta U_i}`

The derivation of the acceptance probability follows that as described in both
TrialComputeMove and ComputeAddAVB with the following transition probabilities.

+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Select a particle of type j and a    |:math:`1/(N_j N_k)`                     |
|particle of type k                   |                                        |
|                                     |                                        |
|[Select a particle of type j and a   |:math:`[1/(N_j N_k)]`                   |
|particle of type k]                  |                                        |
+-------------------------------------+----------------------------------------+
|Select j :math:`\rightarrow` k trial |:math:`1 - P_{bias}`                    |
|type                                 |                                        |
|                                     |                                        |
|[Select k :math:`\rightarrow` j      |:math:`[P_{bias}]`                      |
|trial type]                          |                                        |
+-------------------------------------+----------------------------------------+
|Select site_index_a in AV of         |:math:`1/N^{s,AV,j}_a`                  |
|site_index_j                         |                                        |
|                                     |                                        |
|[Select site_index_a in AV of        |:math:`[1/(N^{s,AV,k}_a + \Delta_{ak})`]|
|site_index_k                         |                                        |
+-------------------------------------+----------------------------------------+
|Place site_index_a in AV of          |:math:`1/v^k_{AV}`                      |
|site_index_k                         |                                        |
|                                     |                                        |
|[Place site_index_a in AV of         |:math:`[1/v^j_{AV}]`                    |
|site_index_j]                        |                                        |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

where :math:`N_t` is the number of particles of type :math:`t`,
:math:`t` is either :math:`j` or :math:`k`,
:math:`v^t_{AV}` is the aggregation volume of particles of type t,
:math:`N^{s,AV,t}_a` is the number of sites with site_index_a, in particles of
type a, that are in the AV of site_index_t,
and :math:`\Delta_{ak}=1`, except when site_index_a is in the AV of site_index_k.
In that special case, :math:`\Delta_{ak}=0`

Application of local detailed balance yields the following acceptance
probability for the :math:`j \rightarrow k` move:

:math:`\chi_ =
\frac{P_{bias}}{1 - P_{bias}}
\frac{N^{s,AV,j}_a}{N^{s,AV,k}_a + \Delta_{ak}}
\frac{v^k_{AV}}{v^j_{AV}}
e^{-\beta\Delta U}`.

The reverse, :math:`k \rightarrow j` acceptance probability is obtained by
simply perturbing the :math:`j,k` indices and inverting the bias probability
ratio as follows:

:math:`\chi_ =
\frac{1 - P_{bias}}{P_{bias}}
\frac{N^{s,AV,k}_a}{N^{s,AV,j}_a + \Delta_{aj}}
\frac{v^j_{AV}}{v^k_{AV}}
e^{-\beta\Delta U}`.

For ease of implementation, we confine ourselves to the case where j and k are
selected from the same pool of particles (e.g., types, or group, or all).
Thus, there is no need to select j->k or k->j.
And because j and k are equivalent types, both terms which involve
:math:`P_{bias}` and terms which involve :math:`v_{AV}` do not contribute.
The resulting acceptance probability for this j->k trial is

:math:`\chi_ =
\frac{N^{s,AV,j}_a}{N^{s,AV,k}_a + \Delta_{ak}}
e^{-\beta\Delta U}`.

\endrst
 */
class TrialAVB4 : public Trial {
 public:
  //@{
  /** @name Arguments
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
    - Trial arguments.
    - TrialStage arguments.
    - SelectParticleAVB arguments.
   */
  explicit TrialAVB4(argtype args = argtype());
  explicit TrialAVB4(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAVB4>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAVB4>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAVB4(std::istream& istr);
  virtual ~TrialAVB4() {}
  //@}
};

inline std::shared_ptr<TrialAVB4> MakeTrialAVB4(argtype args = argtype()) {
  return std::make_shared<TrialAVB4>(args); }

void gen_avb4_args_(argtype * args);

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_AVB4_H_

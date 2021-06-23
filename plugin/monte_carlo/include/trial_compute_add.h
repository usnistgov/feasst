
#ifndef FEASST_MONTE_CARLO_TRIAL_COMPUTE_ADD_H_
#define FEASST_MONTE_CARLO_TRIAL_COMPUTE_ADD_H_

#include <memory>
#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

// HWH add CB+https://hhatch.com/papers/C6SM00473C.pdf
/**
Attempt to add a particle.

The derivation of the acceptance criteria follows a similar procedure as
descibed in TrialComputeMove, except with the following differences.

\rst
The limiting distribution in the grand canonical ensemble is

:math:`\pi_i \propto \frac{e^{-\beta U + \beta \mu_t*N_t}}{\Lambda^{dN}}`

where :math:`\mu_t` is the chemical potential of particles of type t,
:math:`\Lambda` is the de Broglie wavelength, :math:`N_t` is the number
of particles of type t and d is the dimension.

The transition probabilities are as follows, assuming that this move is
coupled with a trial that removes particles with the same selection weight.

+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Select insert trial                  |:math:`1/w`                             |
|                                     |                                        |
|[Select remove trial]                |:math:`[1/w]`                           |
+-------------------------------------+----------------------------------------+
|Place particle of type t             |:math:`1/V`                             |
|                                     |                                        |
|[Delete particle type t]             |:math:`\left[\frac{1}{N_t+1}\right]`    |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

Application of local detailed balance yields the acceptance probability,
:math:`\chi`.

:math:`\frac{e^{-\beta U_o + \beta\mu_t N_t}}{\Lambda^{dN}w V }min(1, \chi) =
\frac{e^{-\beta U_n + \beta\mu_t (N_t+1)}}{\Lambda^{d(N+1)}w (N_t+1)}
min(1, 1/\chi)`

:math:`\chi = \frac{V e^{-\beta\Delta U + \beta\mu_t}}{(N_t+1)\Lambda^d}`

Note that the number of particles, :math:`N_t` is from the perspective of the
old state.
Thus, if the particle has already been added during computation of :math:`\chi`,
then :math:`N_t + 1 \rightarrow N_t`.
The same applies for
:cpp:class:`TrialComputeRemove <feasst::TrialComputeRemove>`.

The de Broglie wavelength, :math:`\Lambda^d`, is absorbed into the definition of
:math:`\mu` for convenience, :math:`\mu + \ln(\Lambda^d)/\beta \rightarrow \mu`.

For configurational bias, consider multiple trial positions and select one.

+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Generate k positions in V.           |:math:`k/V`                             |
|Probability that x_n is in k.        |                                        |
|                                     |                                        |
|[Select particle of type t]          |:math:`[\frac{1}{N_t + 1}]`             |
+-------------------------------------+----------------------------------------+
|Pick x_n in k positions with         |:math:`P_k`                             |
|probability P_k.                     |                                        |
|                                     |                                        |
|[Remove selected particle]           |:math:`[1]`                             |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

:math:`\frac{k P_k}{\Lambda^{dN} V}e^{-\beta U_o + \beta\mu_t N_t} min(1, \chi) =
\frac{1}{\Lambda^{d(N+1)} (N_t+1)}e^{-\beta (U_o + U) + \beta\mu_t (N_t+1)}min(1, 1/\chi)`

where :math:`U` is the interaction energy of the new site with the existing sites and :math:`U_o` is the energy of the original configuration.

:math:`\chi = \frac{V}{k P_k(N_t+1)\Lambda^d}e^{-\beta U + \beta\mu_t}`

If the probability of picking a position is chosen by the Rosenbluth factor,

:math:`P_k = e^{-\beta U}/\sum_i^k e^{-\beta U_i}`.

Thus, the acceptance is given by

:math:`\chi = \frac{V \sum_i^k e^{-\beta U_i}}{k(N_t+1)\Lambda^d}e^{\beta\mu_t}`

For dual-cut configurational bias, the new trials are instead chosen from a
reference potential, :math:`U^r`, that is ideally much faster to compute than
the full potential but still contains sampling-relevant terms (e.g., excluded
volume in a dense system).

:math:`P_k = e^{-\beta U^r}/\sum_i^k e^{-\beta U_i^r}`

:math:`\chi = \frac{V  \sum_i^k e^{-\beta U_i^r}}{k(N_t+1)\Lambda^d}e^{-\beta(U - U^r) + \beta\mu_t}`

Note that these equations consider only a single-site particle.

For the deletion trial, the forward and reverse moves are switched.

+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Select particle of type t.           |:math:`\frac{1}{N_t}`                   |
|                                     |                                        |
|[Generate k positions in V.          |:math:`[k/V]`                           |
|Probability that x_o is in k.]       |                                        |
+-------------------------------------+----------------------------------------+
|Remove selected particle             |:math:`1`                               |
|                                     |                                        |
|[Pick x_o in k positions with        |:math:`[P_k]`                           |
|probability P_k.]                    |                                        |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

:math:`\frac{1}{N_t \Lambda^{dN}}e^{-\beta U_o + \beta\mu_t N_t} min(1, \chi) =
\frac{k P_k}{V \Lambda^{d(N-1)}}e^{-\beta (U_o - U) + \beta\mu_t (N_t-1)}min(1, 1/\chi)`

where :math:`U` is the interaction energy of the removed site with the existing sites and :math:`U_o` is the energy of the original configuration.

:math:`\chi = \frac{k P_k \Lambda^d}{V N_t}e^{\beta U - \beta\mu_t}`

If the probability of picking a position is chosen by the Rosenbluth factor,

:math:`P_k = e^{-\beta U}/\sum_i^k e^{-\beta U_i}`.

Note that one of the k positions is :math:`x_o`.
Thus, the acceptance is given by

:math:`\chi = \frac{k \Lambda^d}{V N_t \sum_i^k e^{-\beta U_i}}e^{-\beta\mu_t}`

\endrst
 */
class TrialComputeAdd : public TrialCompute {
 public:
  explicit TrialComputeAdd(argtype args = argtype());
  explicit TrialComputeAdd(argtype * args);

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeAdd(std::istream& istr);
  virtual ~TrialComputeAdd() {}

 protected:
  void serialize_trial_compute_add_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialComputeAdd> MakeTrialComputeAdd(
    argtype args = argtype()) {
  return std::make_shared<TrialComputeAdd>(args);
}
}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_COMPUTE_ADD_H_

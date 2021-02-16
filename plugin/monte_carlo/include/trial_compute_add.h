
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

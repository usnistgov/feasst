
#ifndef FEASST_MONTE_CARLO_TRIAL_COMPUTE_VOLUME_H_
#define FEASST_MONTE_CARLO_TRIAL_COMPUTE_VOLUME_H_

#include <map>
#include <string>
#include <memory>
#include <vector>
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
Attempt to change the volume

\rst

The derivation of the acceptance criteria follows a similar procedure as
descibed in TrialComputeMove, except with the following differences.
See Lecture 9 of Kofke's CHE 530 course.
Note that, following https://doi.org/10.1080/00268970210126619, this
implementation makes use of a shell particle.
The ideal gas with this algorithm will be :math:`\beta P\langle V\rangle=N`.

The limiting distribution in the isothermal-isobaric ensemble is

:math:`\pi_i \propto e^{-\beta U - \beta P V + (N-1) \ln V} dV`

But if volume is changed in \ln units

:math:`\pi_i \propto e^{-\beta U - \beta P V + N \ln V} d\ln V`

The transition probabilities are as follows, assuming that this move is
coupled with a trial that removes particles with the same selection weight.

+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Select new volume                    |:math:`1/\Delta V`                      |
|                                     |                                        |
|[Select old volume]                  |:math:`[1/\Delta V]`                    |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

Application of local detailed balance yields the acceptance probability,
:math:`\chi`.

:math:`\frac{e^{-\beta U_o - \beta P V_o + (N-1) \ln V_o}}{\Delta V}min(1, \chi) =
\frac{e^{-\beta U_n - \beta P V_n + (N-1) \ln V_n}}{\Delta V}min(1, 1/\chi)`

:math:`\chi = e^{-\beta\Delta U - \beta P \Delta V + (N-1) \ln V_n/V_o}`

And if volume is changed in \ln units

:math:`\chi = e^{-\beta\Delta U - \beta P \Delta V + N \ln V_n/V_o}`
\endrst
 */
class TrialComputeVolume : public TrialCompute {
 public:
  explicit TrialComputeVolume(argtype args = argtype());

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeVolume(std::istream& istr);
  virtual ~TrialComputeVolume() {}

 protected:
  void serialize_trial_compute_add_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialComputeVolume> MakeTrialComputeVolume(
    argtype args = argtype()) {
  return std::make_shared<TrialComputeVolume>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_COMPUTE_VOLUME_H_

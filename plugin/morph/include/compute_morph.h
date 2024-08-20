
#ifndef FEASST_GROWTH_EXPANDED_COMPUTE_MORPH_H_
#define FEASST_GROWTH_EXPANDED_COMPUTE_MORPH_H_

#include <memory>
#include <vector>
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

/**
\rst

Attempt to morph particle(s) a, b, ..., h into particle(s)
:math:`z_a, z_b, ..., z_h`.

To avoid trivial morphs, the particle type of :math:`i` cannot be equal
to the type of :math:`z_i`.

The derivation of the acceptance criteria follows a similar procedure as
descibed in TrialComputeAdd and ComputeAddMultiple.

+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Change particle of type a to z_a     |:math:`1/N_a`                           |
|                                     |                                        |
|[Change particle of type z_a to a]   |:math:`[1/(N_{z_a}+\sum_{x=a}^h         |
|                                     |\delta_{z_az_x} - \delta_{z_ax})]`      |
+-------------------------------------+----------------------------------------+
|Change particle of type b to z_b     |:math:`1/(N_b+1+\sum_{x=a}^{b}          |
|                                     |\delta_{z_x b} - \delta_{xb})`          |
|                                     |                                        |
|[Change particle of type z_b to b]   |:math:`[1/(N_{z_b}+\sum_{x=b}^h         |
|                                     |\delta_{z_bz_x} - \delta_{z_bx})]`      |
+-------------------------------------+----------------------------------------+
|...                                  |...                                     |
+-------------------------------------+----------------------------------------+
|Change particle of type h to z_h     |:math:`1/(N_h+1+\sum_{x=a}^{h}          |
|                                     |\delta_{z_x h} - \delta_{xh})`          |
|                                     |                                        |
|[Change particle of type z_h to h]   |:math:`[1/(N_{z_h} + 1)]`               |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

Application of local detailed balance yields the acceptance probability,

:math:`\chi = e^{-\beta\Delta U} \prod_{i=a}^h
e^{\beta(\mu_{z_i}-\mu_i)}
\frac{
N_i + 1 + \sum_{j=a}^i \delta_{z_j i} - \delta_{ji}
}{
N_{z_i} + \sum_{j=i}^h \delta_{z_i z_j} - \delta_{z_ij}
}`

Note that the number of particles is from the perspective of the old state.

A rederivation using the new state number of particles instead is as follows:

+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Change particle of type a to z_a     |:math:`1/(N_a+\sum_{x=a}^h \delta_{xa}  |
|                                     |- \delta_{z_xa})`                       |
|                                     |                                        |
|[Change particle of type z_a to a]   |:math:`[1/N_{z_a}]`                     |
+-------------------------------------+----------------------------------------+
|Change particle of type b to z_b     |:math:`1/(N_b+\sum_{x=b}^h \delta_{xb}  |
|                                     |- \delta_{z_xb})`                       |
|                                     |                                        |
|[Change particle of type z_b to b]   |:math:`[1/(N_{z_b}+1+\sum_{x=a}^b       |
|                                     |\delta_{xz_b}-\delta_{z_xz_b})]`        |
+-------------------------------------+----------------------------------------+
|...                                  |...                                     |
+-------------------------------------+----------------------------------------+
|Change particle of type h to z_h     |:math:`1/(N_h + 1)`                     |
|                                     |                                        |
|[Change particle of type z_h to h]   |:math:`[1/(N_{z_h}+1+\sum_{x=a}^h       |
|                                     |\delta_{xz_h}-\delta_{z_xz_h})]`        |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

:math:`\chi = e^{-\beta\Delta U} \prod_{i=a}^h
e^{\beta(\mu_{z_i}-\mu_i)}
\frac{
N_i + \sum_{j=i}^h\delta_{ji}-\delta_{z_j i}
}{
N_{z_i}+1+\sum_{j=a}^i\delta_{jz_i}-\delta_{z_jz_i}
}`

\endrst
 */
class ComputeMorph : public TrialCompute {
 public:
  ComputeMorph();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeMorph(std::istream& istr);
  virtual ~ComputeMorph() {}

 protected:
  void serialize_trial_compute_add_(std::ostream& ostr) const;

 private:
  std::vector<int> delta_top_, delta_bottom_;
};

inline std::shared_ptr<ComputeMorph> MakeComputeMorph() {
  return std::make_shared<ComputeMorph>();
}
}  // namespace feasst

#endif  // FEASST_GROWTH_EXPANDED_COMPUTE_MORPH_H_

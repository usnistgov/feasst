
#ifndef FEASST_MAYER_TRIAL_H_
#define FEASST_MAYER_TRIAL_H_

#include "monte_carlo/include/trial_translate.h"

namespace feasst {

class TrialComputeMoveMayer : public TrialCompute {
 public:
  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override {
    compute_rosenbluth(0, criteria, system, acceptance, stages, random);
    acceptance->set_energy_new(acceptance->energy_new());
  }
};

/// Attempt a translation of a random particle, optimized for Mayer-sampling.
class TrialTranslateMayer : public TrialTranslate {
 public:
  TrialTranslateMayer(const argtype& args = argtype()) : TrialTranslate(args) {
    set_mayer(true);
    set(std::make_shared<TrialComputeMoveMayer>());
  }
};

inline std::shared_ptr<TrialTranslateMayer> MakeTrialTranslateMayer(
    const argtype &args = argtype()) {
  return std::make_shared<TrialTranslateMayer>(args);
}

}  // namespace feasst

#endif  // FEASST_MAYER_TRIAL_H_

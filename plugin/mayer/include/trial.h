
#ifndef FEASST_MAYER_TRIAL_H_
#define FEASST_MAYER_TRIAL_H_

#include "monte_carlo/include/trials.h"

namespace feasst {

class TrialComputeMoveNewOnly : public TrialCompute {
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

/// Attempt a translation of a random particle, optimized for not computing old
/// configuration.
inline std::shared_ptr<Trial> MakeTrialTranslateNewOnly(
    const argtype &args = argtype()) {
  auto trial = MakeTrialTranslate(args);
  trial->set_description("TrialTranslateNewOnly");
  trial->set_new_only(true);
  trial->set(std::make_shared<TrialComputeMoveNewOnly>());
  return trial;
}

/// Attempt a rotation of a random particle, optimized for not computing old
/// configuration.
inline std::shared_ptr<Trial> MakeTrialRotateNewOnly(
    const argtype &args = argtype()) {
  auto trial = MakeTrialRotate(args);
  trial->set_description("TrialRotateNewOnly");
  trial->set_new_only(true);
  trial->set(std::make_shared<TrialComputeMoveNewOnly>());
  return trial;
}

}  // namespace feasst

#endif  // FEASST_MAYER_TRIAL_H_


#ifndef FEASST_MAYER_TRIAL_H_
#define FEASST_MAYER_TRIAL_H_

#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"

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
class TrialTranslateNewOnly : public TrialTranslate {
 public:
  TrialTranslateNewOnly(const argtype& args = argtype()) : TrialTranslate(args) {
    set_new_only(true);
    set(std::make_shared<TrialComputeMoveNewOnly>());
  }
};

inline std::shared_ptr<TrialTranslateNewOnly> MakeTrialTranslateNewOnly(
    const argtype &args = argtype()) {
  return std::make_shared<TrialTranslateNewOnly>(args);
}

/// Attempt a rotation of a random particle, optimized for not computing old
/// configuration.
class TrialRotateNewOnly : public TrialRotate {
 public:
  TrialRotateNewOnly(const argtype& args = argtype()) : TrialRotate(args) {
    set_new_only(true);
    set(std::make_shared<TrialComputeMoveNewOnly>());
  }
};

inline std::shared_ptr<TrialRotateNewOnly> MakeTrialRotateNewOnly(
    const argtype &args = argtype()) {
  return std::make_shared<TrialRotateNewOnly>(args);
}

}  // namespace feasst

#endif  // FEASST_MAYER_TRIAL_H_

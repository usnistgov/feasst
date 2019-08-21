
#ifndef FEASST_CHAIN_TRIAL_H_
#define FEASST_CHAIN_TRIAL_H_

#include <vector>
#include <numeric>
#include <string>
#include <memory>
#include "monte_carlo/include/trial.h"
#include "chain/include/trial_select.h"
#include "chain/include/perturb.h"

namespace feasst {

class TrialPivot : public TrialMove {
 public:
  TrialPivot(
    /// These arguments are sent to both PerturbPivot and TrialStage.
    const argtype& args = argtype())
    : TrialMove(
      std::make_shared<TrialSelectEndSegment>(),
      std::make_shared<PerturbPivot>(args),
      args
    ) {
    class_name_ = "TrialPivot";
  }
  virtual ~TrialPivot() {}
};

inline std::shared_ptr<TrialPivot> MakeTrialPivot(
    const argtype &args = argtype()) {
  return std::make_shared<TrialPivot>(args);
}

class TrialCrankshaft : public TrialMove {
 public:
  TrialCrankshaft(
    /// These arguments are sent to both PerturbCrankshaft and TrialStage.
    const argtype& args = argtype())
    : TrialMove(
      std::make_shared<TrialSelectSegment>(),
      std::make_shared<PerturbCrankshaft>(args),
      args
    ) {};

  virtual ~TrialCrankshaft() {}

 protected:
  std::string class_name_ = "TrialCrankshaft";
};

inline std::shared_ptr<TrialCrankshaft> MakeTrialCrankshaft(
    const argtype &args = argtype()) {
  return std::make_shared<TrialCrankshaft>(args);
}

class TrialReptate : public TrialMove {
 public:
  TrialReptate(
    /// These arguments are sent to both PerturbReptate and TrialStage.
    const argtype& args = argtype())
    : TrialMove(
      std::make_shared<TrialSelectReptate>(),
      std::make_shared<PerturbReptate>(args),
      args
    ) {
    class_name_ = "TrialReptate";
  }
  virtual ~TrialReptate() {}
};

inline std::shared_ptr<TrialReptate> MakeTrialReptate(
    const argtype &args = argtype()) {
  return std::make_shared<TrialReptate>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_H_

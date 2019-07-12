
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
    ) {};
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
    ) {};
};

inline std::shared_ptr<TrialReptate> MakeTrialReptate(
    const argtype &args = argtype()) {
  return std::make_shared<TrialReptate>(args);
}

class TrialRegrowLinear : public Trial {
 public:
  TrialRegrowLinear(
    /**
      particle_type : type of particle in configuration (default: 0).
     */
    const argtype& args) : Trial(args) {
    stored_args_ = args;
    set(std::make_shared<TrialComputeMove>());
  }

  void precompute(Criteria * criteria, System * system) override {
    Arguments tmp_args(stored_args_);
    tmp_args.dont_check();
    const int type = tmp_args.key("particle_type").dflt("0").integer();
    const int num_sites = system->configuration().particle_type(type).num_sites();

    // put the first site anywhere
    add_stage(
      std::make_shared<TrialSelectSiteInParticleType>(stored_args_),
      std::make_shared<PerturbAnywhere>(),
      stored_args_
    );

    // for the rest, grow based on bond length only
    for (int site = 1; site < num_sites; ++site) {
      argtype args = stored_args_;
      args.insert(std::pair<std::string, std::string>("mobile_site", str(site)));
      args.insert(std::pair<std::string, std::string>("anchor_site", str(site - 1)));
      add_stage(
        std::make_shared<TrialSelectBond>(args),
        std::make_shared<PerturbDistanceFromAnchor>(args),
        args
      );
    }

    // precompute stages
    Trial::precompute(criteria, system);
  }

 private:
  argtype stored_args_;
};

inline std::shared_ptr<TrialRegrowLinear> MakeTrialRegrowLinear(
    const argtype &args = argtype()) {
  return std::make_shared<TrialRegrowLinear>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_H_

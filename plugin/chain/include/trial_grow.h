
#ifndef FEASST_CHAIN_TRIAL_GROW_H_
#define FEASST_CHAIN_TRIAL_GROW_H_

#include <vector>
#include <string>
#include <memory>
#include "chain/include/trial.h"

namespace feasst {

/**
  Grow a freely jointed linear chain in multiple stages.
 */
class TrialGrowLinear : public Trial {
 public:
  TrialGrowLinear(
    std::shared_ptr<TrialCompute> compute,
    /**
      particle_type : type of particle in configuration (default: 0).
     */
    const argtype& args = argtype()) : Trial(args) {
    stored_args_ = args;
    set(compute);
  }

  void precompute(Criteria * criteria, System * system) override {
    Arguments tmp_args(stored_args_);
    tmp_args.dont_check();
    const int type = tmp_args.key("particle_type").dflt("0").integer();
    const int num_sites = system->configuration().particle_type(type).num_sites();

    // put the first site anywhere
    argtype first_select_args = stored_args_;
    first_select_args.insert({"site", "0"});
    add_stage(
      std::make_shared<TrialSelectParticle>(first_select_args),
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
        std::make_shared<PerturbDistance>(args),
        args
      );
    }

    // precompute stages
    Trial::precompute(criteria, system);
  }

 private:
  argtype stored_args_;
};

inline std::shared_ptr<TrialGrowLinear> MakeTrialGrowLinear(
    std::shared_ptr<TrialCompute> compute,
    const argtype &args = argtype()) {
  return std::make_shared<TrialGrowLinear>(compute, args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_GROW_H_

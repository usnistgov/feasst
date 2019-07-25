
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

// HWH consider using chemical potentials for growth.
// For example. mu = mu_entire_particle/num_growth_stages
// otherwise flat histogram takes care of it just fine
class TrialComputeGrow : public TrialCompute {
 public:
  void perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages) override {
    DEBUG("TrialComputeGrow");
    int old = 0;
//    double factor = 1.;
    if (shrink_ == 1) {
      old = 1;
//      factor = -1.;
    } else if (shrink_ != 0) {
      ERROR("shrink must be initialized: " << shrink_);
    }
    compute_rosenbluth(old, criteria, system, acceptance, stages);
//    const double delta_energy = factor*acceptance->energy_new();
//    DEBUG("current " << criteria->current_energy());
//    DEBUG("delta " << delta_energy);
    acceptance->set_energy_new(criteria->current_energy()
      + acceptance->energy_new()
      - acceptance->energy_old()
    );
//    DEBUG("new " << acceptance->energy_new());
  }

  void set_shrink() { shrink_ = 1; }
  void set_grow() { shrink_ = 0; }

 private:
  int shrink_ = -1;
};

/**
  1-stage growths.
  Each stage attempt makes site physical.
  HWH make a separate trial move that reorders chain sites. This allows either-end trial moves.
    - non-optimal book keeping with cell/neighbor lists.
    - instead, deal with problem like reptation?
  When a growth is attempted with only full chains, attempt creation of a new
    particle and make one site.
  For the reverse, make site unphysical. If only full chains, select
    random particle.
  Precompute all stages, but trial vector has only one and swap them out?
  Macrostate? Use criteria variable "growth_stage?" Or particle itself? Or configuration?
    - selection within configuration, perhaps. try not to add baggage to config
    - criteria "trial_stage" is probably the most general solution.
  Does growth expanded ensemble use chemical potential? Yes, I think it should for each stage.
  Consider using TrialGrow here to manage all the stages and user input of said stages.
  HWH: currently required that first site selection is 0 for shrink
 */
class TrialGrowthExpanded : public Trial {
 public:
  TrialGrowthExpanded(
    std::shared_ptr<Trial> grow,
    // HWH improve interface by using serialization to make a deep copy.
    // then user wouldn't have to provide shrink explicitly
    // as currently implemented, shrink is assumed to be like reverse TrialGrow
    // the selections of first and last stages need to be swapped but accounting
    //  for different sites
    std::shared_ptr<Trial> shrink,
    const argtype& args = argtype())
    : Trial(args) {
    grow_ = grow;
    shrink_ = shrink;
    compute_add_ = std::make_shared<TrialComputeAdd>();
    compute_remove_ = std::make_shared<TrialComputeRemove>();
    compute_grow_ = std::make_shared<TrialComputeGrow>();
    add_(grow_->stages()[0]);
    DEBUG(shrink_);
    DEBUG("num " << shrink_->stages().size());
    const TrialSelect * sel = shrink_->stages()[0]->trial_select();
    // HWH hardcoded for site0
    growing_particle_ = MakeTrialSelectParticle({
      //*(shrink_->stages()[0]->trial_select())
      {"particle_type", str(sel->particle_type())},
      {"site", "0"},
    });

    // the first selection stage of shrink should use perturbed
    auto stage = shrink_->stages()[0];
    stage->set(std::make_shared<TrialSelectPerturbed>());
    shrink_->set(0, stage);
  }

  /// The growth stage relates directly to TrialGrow stages.
  /// Thus, growth stage of 0 refers to the first growth stage where the first
  ///   site of a new particle is being grown.
  /// The last stage is when the particle is completely grown.
  /// However, the reverse move is shifted down by one.
  /// For example, the reverse move at stage 0 would be to select a random
  ///   particle and make the selection of the last stage unphysical
  ///   (e.g., n-1 unphysical)
  /// The reverse move at stage 1 would be to make 0 unphysical.
  /// The reverse move at stage n-1(for n>1) would be to make n-2 unphysical
  int growth_stage() const { return growth_stage_; }

  /// Return the number of stages
  int num_growth_stages() const {
    return static_cast<int>(grow_->num_stages()); }

  // initialize all stage/compute options.
  void precompute(Criteria * criteria, System * system) override {
    grow_->precompute(criteria, system);
    shrink_->precompute(criteria, system);
    // growing_particle_->select(growing_particle_->anchor(), system);
  }

  void before_select(Acceptance * acceptance) override {
    growing_ = random_.coin_flip();
    if (growing_) {
      DEBUG("attempt grow, stage: " << growth_stage_);
      set(0, grow_->stages()[growth_stage_]);
      if (growth_stage_ == 0) {
        set(compute_add_);
      } else {
        compute_grow_->set_grow();
        set(compute_grow_);
      }
    } else {
      DEBUG("attempt shrink, stage: " << growth_stage_);
      if (growth_stage_ == 1) {
        set(compute_remove_);
        set(0, shrink_->stages()[0]);
      } else {
        if (growth_stage_ == 0) {
          set(0, shrink_->stages()[num_growth_stages() - 1]);
        } else {
          set(0, shrink_->stages()[growth_stage_ - 1]);
        }
        compute_grow_->set_shrink();
        set(compute_grow_);
      }
    }

    // acceptanced::perturbed is used for selection of bonds, etc.
    DEBUG("growingp " << growing_particle_->mobile().str());
    acceptance->add_to_perturbed(growing_particle_->mobile());
  }

  // Trial::growth_stage_ -> Criteria::trial_stage?
  bool attempt(Criteria * criteria, System * system) override {
    // whether accepted or rejected, select new growing particle when stage0
    if (growth_stage_ == 0) {
      growing_particle_->select(growing_particle_->anchor(), system);
    }
    const bool accepted = Trial::attempt(criteria, system);
    if (accepted) {
      update_growth_stage_();
      update_growing_particle_();
    }
    if ( (accepted and !growing_) or (!accepted and growing_) ) {
      get_stage_(0)->set_mobile_physical(false, system);
    }
    return accepted;
  }

 private:
  std::shared_ptr<Trial> grow_;
  std::shared_ptr<Trial> shrink_;
  std::shared_ptr<TrialComputeAdd> compute_add_;
  std::shared_ptr<TrialComputeRemove> compute_remove_;
  std::shared_ptr<TrialComputeGrow> compute_grow_;
  std::shared_ptr<TrialSelect> growing_particle_;
  int growth_stage_ = 0;

  // not serialized
  Random random_;
  bool growing_;

  void update_growth_stage_() {
    if (growing_) {
      ++growth_stage_;
      if (growth_stage_ >= num_growth_stages()) {
        growth_stage_ = 0;
      }
    } else {
      --growth_stage_;
      if (growth_stage_ < 0) {
        growth_stage_ = num_growth_stages() - 1;
      }
    }
  }

  void update_growing_particle_() {
    if ( (growing_ and growth_stage_ == 1) or
         (!growing_ and growth_stage_ == num_growth_stages() - 1) ) {
      // select particle that grew/shrunk
      *growing_particle_->get_mobile() = stages()[0]->trial_select()->mobile();

      // HWH hard code growing particle to site 0 for TrialSelectPerturbed
      //  removal. Assumes first site is first stage.
      growing_particle_->get_mobile()->set_site(0, 0, 0);
    }
  }
};

inline std::shared_ptr<TrialGrowthExpanded> MakeTrialGrowthExpanded(
    std::shared_ptr<Trial> grow,
    std::shared_ptr<Trial> shrink,
    const argtype &args = argtype()) {
  return std::make_shared<TrialGrowthExpanded>(grow, shrink, args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_GROW_H_

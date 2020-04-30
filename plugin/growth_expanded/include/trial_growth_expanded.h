
#ifndef FEASST_GROWTH_EXPANDED_TRIAL_GROWTH_EXPANDED_H_
#define FEASST_GROWTH_EXPANDED_TRIAL_GROWTH_EXPANDED_H_

#include <vector>
#include <string>
#include <memory>
#include "chain/include/trial_grow.h"
#include "chain/include/select_perturbed.h"

namespace feasst {

class Random;

class TrialComputeGrowAdd : public TrialCompute {
 public:
  TrialComputeGrowAdd() { class_name_ = "TrialComputeGrowAdd"; }
  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeGrowAdd(std::istream& istr);
  virtual ~TrialComputeGrowAdd() {}
};

class TrialComputeGrowRemove : public TrialCompute {
 public:
  TrialComputeGrowRemove() { class_name_ = "TrialComputeGrowRemove"; }
  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeGrowRemove(std::istream& istr);
  virtual ~TrialComputeGrowRemove() {}
};

// HWH consider using chemical potentials for growth.
// For example. mu = mu_entire_particle/num_growth_stages
// otherwise flat histogram takes care of it just fine
class TrialComputeGrow : public TrialCompute {
 public:
  TrialComputeGrow() { class_name_ = "TrialComputeGrow"; }
  void perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) override {
    DEBUG("TrialComputeGrow");
    ASSERT(shrink_ == 0 or shrink_ == 1,
      "shrink(" << shrink_ << ") must be initialized");
    compute_rosenbluth(shrink_, criteria, system, acceptance, stages, random);
    const TrialSelect * select = (*stages)[0]->trial_select();
    //const int particle_index = select->mobile().particle_index(0);
    //const int particle_type = system->configuration().select_particle(particle_index).type();
    DEBUG("selprob " << select->probability());
    acceptance->set_energy_new(criteria->current_energy()
      + acceptance->energy_new()
      - acceptance->energy_old()
    );

//    const Configuration& config = system->configuration();
//    const double volume = config.domain()->volume();
//    acceptance->add_to_ln_metropolis_prob(
//      log(select->probability())
//      + criteria->beta_mu(particle_type)/criteria->num_trial_states()
//    );
  }

  void set_shrink() { shrink_ = 1; }
  void set_grow() { shrink_ = 0; }

  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeGrow(std::istream& istr);
  virtual ~TrialComputeGrow() {}

 private:
  int shrink_ = -1;
};

/*
Acceptance criteria for growth expanded:
Normally, one has a 1/N probability of selection for delete, and a 1/V prob selection of position
 */

/**
  Separate the growth of a particle into individual Trials.
  Macrostate? Use criteria variable "growth_stage?" Or particle itself? Or configuration?
    - selection within configuration, perhaps. try not to add baggage to config
    - criteria "trial_stage" is probably the most general solution.
  HWH: currently required that first site selection is 0 for shrink
  HWH make a separate trial move that reorders chain sites. This allows either-end trial moves.
    - non-optimal book keeping with cell/neighbor lists.
    - instead, deal with problem like reptation?
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
    class_name_ = "TrialGrowthExpanded";
    grow_ = grow;
    shrink_ = shrink;
    compute_add_ = std::make_shared<TrialComputeGrowAdd>();
    compute_remove_ = std::make_shared<TrialComputeGrowRemove>();
    compute_grow_ = std::make_shared<TrialComputeGrow>();
    add_(grow_->stages()[0]);
    DEBUG(shrink_);
    DEBUG("num " << shrink_->stages().size());
    const TrialSelect * sel = shrink_->stages()[0]->trial_select();
    growing_particle_ = MakeTrialSelectParticle({
      {"particle_type", str(sel->particle_type())},
      {"site", "0"}, // HWH hardcoded for site0
    });

    // the first selection stage of shrink should use perturbed
    auto stage = shrink_->stages()[0];
    stage->set(std::make_shared<SelectPerturbed>());
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

  void precompute(Criteria * criteria, System * system) override {
    grow_->precompute(criteria, system);
    shrink_->precompute(criteria, system);
    criteria->set_trial_state(growth_stage_, num_growth_stages());
  }

  void before_select(Acceptance * acceptance, Criteria * criteria) override {
    criteria->set_trial_state(current_growth_stage_(growing_),
                              num_growth_stages());
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
      if (growth_stage_ == 0) {
        set(0, shrink_->stages()[num_growth_stages() - 1]);
        compute_grow_->set_shrink();
        set(compute_grow_);
      } else if (growth_stage_ == 1) {
        set(0, shrink_->stages()[0]);
        set(compute_remove_);
      } else {
        set(0, shrink_->stages()[growth_stage_ - 1]);
        compute_grow_->set_shrink();
        set(compute_grow_);
      }
    }

    // acceptanced::perturbed is used for selection of bonds, etc.
    DEBUG("growingp " << growing_particle_->mobile().str());
    if (growth_stage_ != 0) {
      acceptance->add_to_perturbed(growing_particle_->mobile());
    }
  }

  // Trial::growth_stage_ -> Criteria::trial_stage?
  bool attempt(Criteria * criteria, System * system, Random * random) override;

  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialGrowthExpanded(std::istream& istr);
  virtual ~TrialGrowthExpanded() {}

 private:
  std::shared_ptr<Trial> grow_;
  std::shared_ptr<Trial> shrink_;
  std::shared_ptr<TrialComputeGrow> compute_grow_;
  std::shared_ptr<TrialComputeGrowAdd> compute_add_;
  std::shared_ptr<TrialComputeGrowRemove> compute_remove_;
  std::shared_ptr<TrialSelect> growing_particle_;
  int growth_stage_ = 0;

  // not serialized
  bool growing_;

  int current_growth_stage_(const bool growing) const {
    int stage = growth_stage_;
    if (growing) {
      ++stage;
      if (stage >= num_growth_stages()) {
        stage = 0;
      }
    } else {
      --stage;
      if (stage < 0) {
        stage = num_growth_stages() - 1;
      }
    }
    return stage;
  }

  void update_growing_particle_() {
    if ( (growing_ and growth_stage_ == 1) or
         (!growing_ and growth_stage_ == num_growth_stages() - 1) ) {
      *growing_particle_->get_mobile() = stages()[0]->trial_select()->mobile();
      // HWH hard code growing particle to site 0 for SelectPerturbed
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

#endif  // FEASST_GROWTH_EXPANDED_TRIAL_GROWTH_EXPANDED_H_


#ifndef FEASST_CHAIN_TRIAL_REGROW_H_
#define FEASST_CHAIN_TRIAL_REGROW_H_

#include "monte_carlo/include/trial_translate.h"
#include "chain/include/perturb_regrow.h"

namespace feasst {

/**
  Assumes a linear chain.
 */
class StageFactoryRegrow : public StageFactory {
 public:
  void parse_select(const SelectList& select) override {
    set_max_stage(select.num_sites() - 1);

    DEBUG("num stages " << num_stages());
    bool reverse = false;
    DEBUG("select " << select.str());
    if (select.site_indices()[0][0] == 0) {
      reverse = true;
    }

    SelectList rebond = select, anchor = rebond;
    if (reverse) {
      rebond.remove_last_site();
      anchor.remove_first_site();
    } else {
      rebond.remove_first_site();
      anchor.remove_last_site();
    }

    Select excluded = rebond;

    for (int stage_index = 0; stage_index < rebond.num_sites(); ++stage_index) {
      SelectList parsed_rb = rebond;
      SelectList parsed_an = anchor;
      int pop_fronts = stage_index;
      int pop_backs = rebond.num_sites() - 1 - stage_index;
      int index = pop_fronts;
      if (reverse) {
        swap(&pop_fronts, &pop_backs);
        index = pop_backs;
      }
      DEBUG("pop_backs/fronts " << pop_backs << " " << pop_fronts);
      ASSERT(pop_fronts + pop_backs == rebond.num_sites() - 1,
        "linear regrow stages by 1 sites each.");
      for (int pop_back = 0; pop_back < pop_backs; ++pop_back) {
        parsed_rb.remove_last_site();
        parsed_an.remove_last_site();
      }
      for (int pop_front = 0; pop_front < pop_fronts; ++pop_front) {
        parsed_rb.remove_first_site();
        parsed_an.remove_first_site();
      }
      DEBUG("pop_backs " << pop_backs);
      DEBUG("in parsed " << parsed_rb.num_sites());
      DEBUG("index " << index);

      // add the excluded sites (e.g., not yet regrown)
      excluded.remove(parsed_rb);
      DEBUG("to be excluded " << excluded.str());
      parsed_rb.exclude(excluded);

      stages_[index]->parse_select(parsed_rb);
      stages_[index]->parse_anchor(parsed_an);
    }
  }

  void update_select(std::shared_ptr<Stage> stage, const System * system) override {
    SelectList select = stage->perturb()->selection();
    SelectList anchor = stage->perturb()->anchor();
    DEBUG("moving " << select.str());
    select.load_positions(system->configuration().particles());
    anchor.load_positions(system->configuration().particles());
    stage->parse_select(select);
    stage->parse_anchor(anchor);
  }
};

/**
  For a staged regrow,
  - variable size stages? resize if larger requested. Truncate to max stage for given trial?
  - each stage is has perturb_regrow, where selection is parsed for just two sites: a dangling site and the fixed one its bonded to
  - assume linear chain to set up perturb selection parsing
  - intra interactions must ignore those sites which have yet to be regrown
 */
class TrialRegrow : public TrialMove, public StagedTrial {
 public:
  TrialRegrow(
    /**
      reference : index of the reference potential.
      num_steps : number of steps per stage.
      max_length : maximum length of selected segment. If -1 (default), then
        randomly select all possible lengths.
     */
    const argtype &args = argtype()) : TrialMove(args) {
    regrow_ = std::make_shared<PerturbRegrow>();
    set_perturb(regrow_);
    DEBUG("reg " << regrow_.get()); //
    args_.init(args);
    max_length_ = args_.key("max_length").dflt("-1").integer();
    auto stage = std::make_shared<Stage>();
    parse_ref_and_num_steps(args, &args_, stage);
    stage->set(std::make_shared<PerturbRegrow>(*regrow_));
    stages_.add(stage);
  }

  /// Initialize the number of stages to the number of site types.
  void precompute(const std::shared_ptr<Criteria> criteria,
    const System& system) override {
    /// find the particle with the largest number of sites
    int max = 0;
    for (const Particle part : system.configuration().particle_types().particles()) {
      if (part.num_sites() > max) {
        max = part.num_sites();
      }
    }
    DEBUG("max sites " << max);

    /// set this many stages, minus the already existing one.
    for (int type = 0; type < max - 1; ++type) {
      auto stg = std::make_shared<Stage>(*stages_.stage(0));
      stg->set(std::make_shared<PerturbRegrow>(*regrow_));
      stages_.add(stg);
    }
  }

  void select(System * system) override {
    regrow_->select_random_end_segment_in_particle(
      group_index(),
      system->configuration(),
      max_length_);
    DEBUG("selection to regrow: " << regrow_->selection().str());
    if (regrow_->selection().num_sites() != 0) {
      stages_.parse_select(regrow_->selection());
    }
  }

  void move_and_acceptance(Criteria * criteria, System * system, AcceptanceCriteria * accept_criteria_) override {
    DEBUG("energy of select " << system->energy(regrow_->selection()));
    DEBUG("here i am");
    regrow_->get_config_before_move(system);  // required for reverting
    int reject = 0;
    double en_old;
    const double ln_met_old = stages_.compute_rosenbluth(
      1,  // old
      regrow_->selection(),
      criteria,
      system,
      &en_old,
      &reject);
    double en_new;
    const double ln_met_new = stages_.compute_rosenbluth(
      0,  // new
      regrow_->selection(),
      criteria,
      system,
      &en_new,
      &reject);
    const double delta_energy = en_new - en_old;
    DEBUG("here i am2");
    DEBUG("reject? " << reject);
    accept_criteria_->force_rejection = reject;
    accept_criteria_->ln_metropolis_prob -= ln_met_old;
    accept_criteria_->ln_metropolis_prob += ln_met_new;
    accept_criteria_->energy_new = criteria->current_energy() + delta_energy;
    accept_criteria_->energy_new_select = en_new;
    accept_criteria_->system = system;
    DEBUG("new en " << accept_criteria_->energy_new);
    DEBUG("delta_energy " << delta_energy);
    DEBUG("ln_met " << accept_criteria_->ln_metropolis_prob);
    DEBUG("sel " << regrow_->selection().str() << " reg " << regrow_.get());
    DEBUG("perturb addresss " << perturb().get());
    DEBUG("perturb sel " << perturb()->selection().str());
    regrow_->after_move();  // required for reverting
  }

 private:
  StageFactoryRegrow stages_;
  std::shared_ptr<PerturbRegrow> regrow_;
  int max_length_;
};

inline std::shared_ptr<TrialRegrow> MakeTrialRegrow(
    const argtype &args = argtype()) {
  return std::make_shared<TrialRegrow>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_REGROW_H_

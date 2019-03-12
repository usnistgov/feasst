
#ifndef FEASST_CORE_ROSENBLUTH_H_
#define FEASST_CORE_ROSENBLUTH_H_

#include "core/include/trial.h"

namespace feasst {

/**
  Multistage:
  -take the original selection: divide it into an ordered list of selections
   each of these are referred to as stages
  -for each of these stages, attempt "n" perturbations on select, storing energy.
  -use energy to select one of the perturbations for a given stage.

  -after all stages, accept or reject using rosenbluth factors, plus other energy terms
  - that may not have been considered (e.g., ref potential?)

  Need
    - turn off interactions for selections which have not been perturbed yet
      - only an issue for the intra-particle?
    - for reference potential use, also may provide a group index?
 */

// perform trial n times using perturb, accept or reject based on rosenbluth.
// optionally, for old configuration, simply revert? (or handle this in higher class)
// ^-- optimization

class Rosenbluth {
 public:
  void resize(const int num) {
    energy_.resize(num);
    boltzman_.resize(num);
    stored_.resize(num);
  }
  int num() const { return static_cast<int>(energy_.size()); }
  void set_energy(const int index, const double energy) {
    energy_[index] = energy;
  }
  void compute(const double beta) {
    for (int index = 0; index < num(); ++index) {
      boltzman_[index] = exp(-beta*energy_[index]);
    }
    TRACE("boltzman " << feasst_str(boltzman_));
    rosenbluth_ = std::accumulate(boltzman_.begin(), boltzman_.end(), 0.);
    TRACE("rosen " << rosenbluth_);
    DEBUG("energy " << feasst_str(energy_));
    if (rosenbluth_ <= 0) {
      index_ = -1;
      return;
    }
    cumulative_ = cumulative_probability(boltzman_);
    TRACE("cumulative " << feasst_str(cumulative_));
    index_ = random_.index_from_cumulative_probability(cumulative_);
  }

  void store(const int index, const SelectList& select, const System * system) {
    stored_[index].store(select, system->configuration());
  }

  const SelectList& chosen() const {
    ASSERT(index_ != -1, "error");
    return stored_[index_];
  }

  double chosen_energy() const {
    DEBUG("chosen index " << index_);
    if (index_ == -1) {
      return energy_[0];
    }
    return energy_[index_];
  }
  double rosenbluth() const { return rosenbluth_; }

  double energy(const int index) const { return energy_[index]; }

  int chosen_index() const { return index_; }

 private:
  std::vector<double> energy_;
  std::vector<double> boltzman_;
  std::vector<double> cumulative_;
  std::vector<SelectList> stored_;
  double rosenbluth_;
  int index_;
  Random random_;
};

class Stage {
 public:
  void parse_select(const SelectList& select) {
    DEBUG("parsing select " << select.str());
    perturb_->set_selection(select);
  }

  void parse_anchor(const SelectList& anchor) {
    perturb_->set_anchor(anchor);
  }

  void compute(const int old, Criteria* criteria, System * system) {
    for (int step = 0; step < rosenbluth_.num(); ++step) {
      // skip first perturbation if old
      if (step != 0 || old != 1) {
        perturb_->perturb(system);
      }
      rosenbluth_.store(step, perturb_->selection(), system);
      if (reference_ == -1) {
        DEBUG("select " << perturb_->selection().str());
        rosenbluth_.set_energy(step, system->energy(perturb_->selection()));
      } else {
        rosenbluth_.set_energy(step,
          system->reference_energy(perturb_->selection(), reference_));
      }
      if (step != 0 || old != 1) {
        perturb_->revert();
      }
    }
    rosenbluth_.compute(criteria->beta());
    if (old != 1 and rosenbluth_.chosen_index() != -1) {
      DEBUG("updating positions " << rosenbluth_.chosen().str());
      DEBUG("pos0 " << rosenbluth_.chosen().site_positions()[0][0].str());
      // DEBUG("pos1 " << rosenbluth_.chosen().site_positions()[0][1].str());
      system->get_configuration()->update_positions(rosenbluth_.chosen());
    }
  }

  void set_reference(const int ref = -1) { reference_ = ref; }
  void set_num_steps(const int num) { rosenbluth_.resize(num); }
  void set(std::shared_ptr<Perturb> perturb) { perturb_ = perturb; }
  const std::shared_ptr<Perturb> perturb() const { return perturb_; }
  const Rosenbluth& rosenbluth() const { return rosenbluth_; }
  int reference() const { return reference_; }

//  Stage deep_copy() const {
//    Stage stage = *this;
//    stage.perturb_ = std::make_shared<Perturb>(*perturb_);
//    return stage;
//  }

 private:
  std::shared_ptr<Perturb> perturb_;
  // selection indices for when taking piece of whole selection
  // SelectList select_;
  int reference_ = -1; // if >= 0, use ref
  Rosenbluth rosenbluth_;
};

// divide up selection
// store/perform each stage
// HWH site: frenkel and smit, dual-cut CB, etc
class StageFactory {
 public:
  StageFactory() {
    set_max_stage();
  }

  /// Add another stage.
  void add(std::shared_ptr<Stage> stage) {
    stages_.push_back(stage);
    update_ref_flag_();
  }

  /// Create another stage like the last one.
  void push_back() {
    stages_.push_back(stages_.back()); }
  //stages_.push_back(std::make_shared<Stage>(stages_.back()->deep_copy())); }

  virtual void parse_select(const SelectList& select) {
    for (std::shared_ptr<Stage> stage : stages_) {
      stage->parse_select(select);
    }
  }

  double compute_rosenbluth(
      const int old,
      const SelectList& select,
      Criteria* criteria,
      System * system,
      double * energy_change,
      int * reject) {
    double ln_rosenbluth = 0.;
    *energy_change = 0.;
    double ln_metropolis_prob = 0.;
    auto stages = stages_;
    if (max_stage_ != -1) {
      stages.resize(max_stage_);
    }
    DEBUG("max stages " << max_stage_ << " size " << stages.size());
    for (std::shared_ptr<Stage> stage : stages) {
      if (old != 1) {
        update_select(stage, system);
      }
      stage->compute(old, criteria, system);
      if (stage->rosenbluth().chosen_index() == -1) {
        *energy_change = NEAR_INFINITY;
        *reject = 1;
        DEBUG("auto reject");
        return NEAR_INFINITY;
      }
      ln_rosenbluth += log(stage->rosenbluth().rosenbluth());
      if (old == 1 && stage->rosenbluth().chosen_index() != -1) {
        *energy_change += stage->rosenbluth().energy(0);
      } else {
        *energy_change += stage->rosenbluth().chosen_energy();
      }
    }
    if (is_ref_used_) {
      const double en_full = system->energy(select);
      ln_metropolis_prob += -1.*criteria->beta()*(en_full - *energy_change);
      *energy_change = en_full;
    }
    ln_metropolis_prob += ln_rosenbluth;
    return ln_metropolis_prob;
  }

  // consider using this in base class instead of parse_select
  virtual void update_select(std::shared_ptr<Stage> stage, const System * system) {}

  /// Limit the computation to a number of stages.
  void set_max_stage(
    /// The default value of -1 for max considers all stages.
    const int max = -1) { max_stage_ = max;}

  int num_stages() const { return static_cast<int>(stages_.size()); }

  const std::shared_ptr<Stage> stage(const int index) const { return stages_[index]; }

 protected:
  std::vector<std::shared_ptr<Stage> > stages_;

 private:
  bool is_ref_used_;
  int max_stage_;

  void update_ref_flag_() {
    is_ref_used_ = false;
    for (const std::shared_ptr<Stage> stage : stages_) {
      if (stage->reference() != -1) {
        is_ref_used_ = true;
      }
    }
  }
};

}  // namespace feasst

#endif  // FEASST_CORE_ROSENBLUTH_H_


#ifndef FEASST_CORE_MODIFY_H_
#define FEASST_CORE_MODIFY_H_

#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include "core/include/analyze.h"

namespace feasst {

class Modify : public Stepper {
 public:
  virtual void initialize(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) {
    // do nothing by default
  }

  virtual void trial(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) {
    if (is_time(steps_per_update(), get_steps_since_update())) {
      update(criteria, system, trial_factory);
    }
    if (is_time(steps_per_write(), get_steps_since_write())) {
      printer(write(criteria, system, trial_factory));
    }
  }

  virtual void update(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) {
    ERROR("not implemented");
  }

  virtual std::string write(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) {
    ERROR("not implemented");
    return std::string("");
  }

  virtual ~Modify() {}
};

class ModifyFactory : public Modify {
 public:
  ModifyFactory() : Modify() {}

  void add(std::shared_ptr<Modify> modify) {
    modifiers_.push_back(modify);
  }

  const std::vector<std::shared_ptr<Modify> >& modifiers() const {
    return modifiers_; }

  void trial(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) override {
    for (const std::shared_ptr<Modify> modify : modifiers_) {
      modify->trial(criteria, system, trial_factory);
    }
  }

  private:
    std::vector<std::shared_ptr<Modify> > modifiers_;
};

class ModifyUpdateOnly : public Modify {
 public:
  ModifyUpdateOnly() : Modify() {
    // disable write
    Modify::set_steps_per_write(-1); }

  void set_steps_per_write(const int steps) override {
    ERROR("This modify is update only."); }

  void set_steps_per(const int steps) { set_steps_per_update(steps); }
};

/**
  Check that the running energy from criteria is equivalent, within tolerance,
  to a fresh (unoptimized) calculation over the entire configuration.
 */
class EnergyCheck : public ModifyUpdateOnly {
 public:
  void set_tolerance(const double tolerance) { tolerance_ = tolerance; }

  void update(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) override {
    const double energy = system->unoptimized_energy();
    const double running_energy = criteria->running_energy();
    ASSERT(std::abs(energy - running_energy) < tolerance_,
      "Energy check failure. There is a problem with the potentials. " <<
      "The unoptimized energy of the entire configuration was computed as " <<
      energy << " but the running energy from criteria (the accumulation of a "
      << "change in energy over a series of steps) is " << running_energy <<
      ". The difference(" << std::abs(energy - running_energy) << ") is " <<
      "greater than the tolerance(" << tolerance_ << ")");
    criteria->set_running_energy(energy);
  }

 private:
  double tolerance_;
};

/**
  Check that the running energy from criteria is equivalent, within tolerance,
  to a fresh (unoptimized) calculation over the entire configuration.
 */
class Tuner : public ModifyUpdateOnly {
 public:
  void update(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) override {
    trial_factory->tune();
  }
};

}  // namespace feasst

#endif  // FEASST_CORE_MODIFY_H_

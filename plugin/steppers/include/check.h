
#ifndef FEASST_STEPPERS_CHECK_H_
#define FEASST_STEPPERS_CHECK_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Run periodic checks on each class.
 */
class Check : public ModifyUpdateOnly {
 public:
  Check(const argtype &args = argtype()) : ModifyUpdateOnly(args) {}
  void update(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) override {
    system->configuration().check();
  }

  std::string class_name() const override { return std::string("Check"); }

  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(1, ostr);
  }

  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<Check>(istr); }

  Check(std::istream& istr) : ModifyUpdateOnly(istr) {
    feasst_deserialize_version(istr); }
};

inline std::shared_ptr<Check> MakeCheck(const argtype &args = argtype()) {
  return std::make_shared<Check>(args);
}

/**
  Check that the running energy from criteria is equivalent, within tolerance,
  to an (unoptimized) calculation over the entire configuration.
 */
class EnergyCheck : public ModifyUpdateOnly {
 public:
  EnergyCheck(
    /**
      tolerance : relative absolute difference between running energy
        and recomputed energy.
    */
    const argtype &args = argtype()) : ModifyUpdateOnly(args) {
    // default
    set_tolerance();

    // parse
    args_.init(args);
    if (!args_.key("tolerance").empty()) {
      set_tolerance(args_.dble());
    }
    check_ = MakeCheck();
  }

  /// Set tolerance for energy drift per check.
  void set_tolerance(const double tolerance = 1e-10) { tolerance_ = tolerance; }

  void update(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) override {
    const double energy = system->unoptimized_energy();
    const double current_energy = criteria->current_energy();
    ASSERT(std::abs(energy - current_energy) < tolerance_,
      "Energy check failure. There is a problem with the potentials. " <<
      "The unoptimized energy of the entire configuration was computed as " <<
      energy << " but the running energy from criteria (the accumulation of a "
      << "change in energy over a series of steps) is " << current_energy <<
      ". The difference(" << std::abs(energy - current_energy) << ") is " <<
      "greater than the tolerance(" << tolerance_ << ")");
    criteria->set_current_energy(energy);
    check_->update(criteria, system, trial_factory);
  }

  std::string class_name() const override { return std::string("EnergyCheck"); }

  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(1, ostr);
    feasst_serialize(tolerance_, ostr);
    feasst_serialize_fstdr(check_, ostr);
  }

  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<EnergyCheck>(istr); }

  EnergyCheck(std::istream& istr) : ModifyUpdateOnly(istr) {
    feasst_deserialize_version(istr);
    feasst_deserialize(&tolerance_, istr);
    // feasst_deserialize_fstdr(modify->check_, istr);
    { // HWH for unknown reasons the above template function does not work
      int existing;
      istr >> existing;
      if (existing != 0) {
        check_ = check_->deserialize(istr);
      }
    }
  }

 private:
  double tolerance_;
  std::shared_ptr<Modify> check_;
};

inline std::shared_ptr<EnergyCheck> MakeEnergyCheck(const argtype &args = argtype()) {
  return std::make_shared<EnergyCheck>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CHECK_H_

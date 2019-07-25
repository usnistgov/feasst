
#ifndef FEASST_STEPPERS_CHECK_ENERGY_H_
#define FEASST_STEPPERS_CHECK_ENERGY_H_

#include "monte_carlo/include/modify.h"
#include "steppers/include/check.h"

namespace feasst {

/**
  Check that the running energy from criteria is equivalent, within tolerance,
  to an (unoptimized) calculation over the entire configuration.
 */
class CheckEnergy : public ModifyUpdateOnly {
 public:
  CheckEnergy(
    /**
      tolerance : relative absolute difference between running energy
        and recomputed energy (default: 1e-10).
    */
    const argtype &args = argtype()) : ModifyUpdateOnly(args) {

    // parse
    args_.init(args);
    tolerance_ = args_.key("tolerance").dflt(str(1e-10)).dble();
    check_ = MakeCheck();
  }

  void update(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override {
    check_->update(criteria, system, trial_factory);
    DEBUG("computing unoptimized energy for check");
    const double energy = system->unoptimized_energy();
    const double current_energy = criteria->current_energy();
    DEBUG("energy " << energy);
    ASSERT(std::abs(energy - current_energy) < tolerance_,
      "Energy check failure. There is a problem with the potentials. " <<
      "The unoptimized energy of the entire configuration was computed as " <<
      energy << " but the running energy from criteria (the accumulation of a "
      << "change in energy over a series of steps) is " << current_energy <<
      ". The difference(" << std::abs(energy - current_energy) << ") is " <<
      "greater than the tolerance(" << tolerance_ << "). "
      << system->unoptimized().str());
    criteria->set_current_energy(energy);
  }

  std::string class_name() const override { return std::string("CheckEnergy"); }

  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(1, ostr);
    feasst_serialize(tolerance_, ostr);
    feasst_serialize_fstdr(check_, ostr);
  }

  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<CheckEnergy>(istr); }

  CheckEnergy(std::istream& istr) : ModifyUpdateOnly(istr) {
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

inline std::shared_ptr<CheckEnergy> MakeCheckEnergy(const argtype &args = argtype()) {
  return std::make_shared<CheckEnergy>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CHECK_ENERGY_H_

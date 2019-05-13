
#ifndef FEASST_STEPPERS_ENERGY_H_
#define FEASST_STEPPERS_ENERGY_H_

#include "monte_carlo/include/analyze.h"
#include "math/include/accumulator.h"

namespace feasst {

/**
  Accumulate average energy.
 */
class Energy : public Analyze {
 public:
  Energy(const argtype &args = argtype());

  void update(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    energy_.accumulate(criteria->current_energy());
  }

  std::string write(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    std::stringstream ss;
    ss << energy_.average() << " ";
    // INFO(ss.str());
    return ss.str();
  }

  std::string class_name() const override { return std::string("Energy"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Energy>(istr); }
  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(1, ostr);
    feasst_serialize_fstobj(energy_, ostr);
  }
  Energy(std::istream& istr) : Analyze(istr) {
    feasst_deserialize_version(istr);
    feasst_deserialize_fstobj(&energy_, istr);
  }

 private:
  Accumulator energy_;
};

inline std::shared_ptr<Energy> MakeEnergy(const argtype &args = argtype()) {
  return std::make_shared<Energy>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_ENERGY_H_

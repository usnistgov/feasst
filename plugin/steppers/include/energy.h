
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
  Energy(
    /**
      num_block : number of updated per block (default: 1e5).
     */
    const argtype &args = argtype());

  void update(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    energy_.accumulate(criteria->current_energy());
  }

  std::string write(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    std::stringstream ss;
    ss << energy_.str() << " ";
    DEBUG(ss.str());
    return ss.str();
  }

  const Accumulator& energy() const { return energy_; }

  const Accumulator& accumulator() const override { return energy(); }

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

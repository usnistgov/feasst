
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
  /**
    args:
    - num_block: number of updated per block (default: 1e5).
   */
  explicit Energy(const argtype &args = argtype());

  void update(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  std::string write(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  const Accumulator& energy() const { return energy_; }

  const Accumulator& accumulator() const override { return energy(); }

  // serialize
  std::string class_name() const override { return std::string("Energy"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Energy>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit Energy(std::istream& istr);

 private:
  Accumulator energy_;
};

inline std::shared_ptr<Energy> MakeEnergy(const argtype &args = argtype()) {
  return std::make_shared<Energy>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_ENERGY_H_

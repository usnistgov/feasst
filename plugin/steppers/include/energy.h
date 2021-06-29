
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
  explicit Energy(argtype args = argtype());
  explicit Energy(argtype * args);

  std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const override;

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  /// Return the energy.
  const Accumulator& energy() const { return accumulator(); }

  // serialize
  std::string class_name() const override { return std::string("Energy"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Energy>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<Energy>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Energy(std::istream& istr);
  explicit Energy(const Analyze& energy);
};

inline std::shared_ptr<Energy> MakeEnergy(argtype args = argtype()) {
  return std::make_shared<Energy>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_ENERGY_H_

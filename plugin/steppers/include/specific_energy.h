
#ifndef FEASST_STEPPERS_SPECIFIC_ENERGY_H_
#define FEASST_STEPPERS_SPECIFIC_ENERGY_H_

#include "monte_carlo/include/analyze.h"
#include "math/include/accumulator.h"

namespace feasst {

/**
  Accumulate average energy per number of particles.
 */
class SpecificEnergy : public Analyze {
 public:
  //@{
  /** @name Arguments
    - Stepper arguments.
   */
  explicit SpecificEnergy(argtype args = argtype());
  explicit SpecificEnergy(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(const MonteCarlo& mc) override;
  std::string write(const MonteCarlo& mc) override;

  /// Return the energy.
  const Accumulator& energy() const { return accumulator(); }

  // serialize
  std::string class_name() const override { return std::string("SpecificEnergy"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<SpecificEnergy>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<SpecificEnergy>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit SpecificEnergy(std::istream& istr);
  //@}
};

}  // namespace feasst

#endif  // FEASST_STEPPERS_SPECIFIC_ENERGY_H_

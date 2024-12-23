
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
  //@{
  /** @name Arguments
    - Stepper arguments.
   */
  explicit Energy(argtype args = argtype());
  explicit Energy(argtype * args);

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
  std::string class_name() const override { return std::string("Energy"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Energy>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<Energy>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Energy(std::istream& istr);
  explicit Energy(const Analyze& energy);
  //@}
};

inline std::shared_ptr<Energy> MakeEnergy(argtype args = argtype()) {
  return std::make_shared<Energy>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_ENERGY_H_


#ifndef FEASST_STEPPERS_HEAT_CAPACITY_H_
#define FEASST_STEPPERS_HEAT_CAPACITY_H_

#include "monte_carlo/include/analyze.h"
#include "math/include/accumulator.h"

namespace feasst {

/**
  Accumulate average heat capacity, C, per Boltzmann constant (\f$k_B\f$),
  as computed from the fluctuation in the energy, U,

  \f$ C/k_B = \beta^2 (\langle U^2 \rangle - \langle U \rangle^2) \f$.

  Note that the heat capacity can also be computed from the moments output by
  Energy.
 */
class HeatCapacity : public Analyze {
 public:
  //@{
  /** @name Arguments
    - Stepper arguments.
   */
  explicit HeatCapacity(argtype args = argtype());
  explicit HeatCapacity(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(const MonteCarlo& mc) override;
  std::string write(const MonteCarlo& mc) override;

  // serialize
  std::string class_name() const override { return std::string("HeatCapacity"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<HeatCapacity>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<HeatCapacity>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit HeatCapacity(std::istream& istr);
  explicit HeatCapacity(const Analyze& energy);
  //@}

 private:
  Accumulator energy_;
};

inline std::shared_ptr<HeatCapacity> MakeHeatCapacity(argtype args = argtype()) {
  return std::make_shared<HeatCapacity>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_HEAT_CAPACITY_H_

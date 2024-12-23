
#ifndef FEASST_STEPPERS_VOLUME_H_
#define FEASST_STEPPERS_VOLUME_H_

#include "monte_carlo/include/analyze.h"
#include "math/include/accumulator.h"

namespace feasst {

/**
  Accumulate average number density.
 */
class Density : public Analyze {
 public:
  //@{
  /** @name Arguments
    - Stepper arguments.
   */
  explicit Density(argtype args = argtype());
  explicit Density(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(const MonteCarlo& mc) override;
  std::string write(const MonteCarlo& mc) override;

  const Accumulator& density() const { return accumulator(); }

  // serialize
  std::string class_name() const override { return std::string("Density"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Density>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<Density>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Density(std::istream& istr);
  //@}
};

inline std::shared_ptr<Density> MakeDensity(argtype args = argtype()) {
  return std::make_shared<Density>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_VOLUME_H_

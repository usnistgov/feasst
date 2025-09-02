
#ifndef FEASST_STEPPERS_SPECIFIC_VOLUME_H_
#define FEASST_STEPPERS_SPECIFIC_VOLUME_H_

#include "monte_carlo/include/analyze.h"
#include "math/include/accumulator.h"

namespace feasst {

/**
  Accumulate average volume per number of particles.
 */
class SpecificVolume : public Analyze {
 public:
  //@{
  /** @name Arguments
    - Stepper arguments.
   */
  explicit SpecificVolume(argtype args = argtype());
  explicit SpecificVolume(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(const MonteCarlo& mc) override;
  std::string write(const MonteCarlo& mc) override;

  const Accumulator& volume() const { return accumulator(); }

  // serialize
  std::string class_name() const override { return std::string("SpecificVolume"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<SpecificVolume>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<SpecificVolume>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit SpecificVolume(std::istream& istr);
  //@}
};

}  // namespace feasst

#endif  // FEASST_STEPPERS_SPECIFIC_VOLUME_H_

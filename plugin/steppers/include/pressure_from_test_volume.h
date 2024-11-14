
#ifndef FEASST_STEPPERS_PRESSURE_FROM_TEST_VOLUME_H_
#define FEASST_STEPPERS_PRESSURE_FROM_TEST_VOLUME_H_

#include <vector>
#include <memory>
#include "monte_carlo/include/modify.h"

namespace feasst {

class Accumulator;

typedef std::map<std::string, std::string> argtype;

/**
  Compute pressure using test volume changes.
  See https://doi.org/10.1063/1.472721

  Hard systems must consider only volume reductions
  (i.e., negative delta_volume).

  \f$p = \frac{1}{\beta \Delta V}\ln \left\langle \left(\frac{V'}{V}\right)^N e^{-\beta\Delta U}\right\rangle \f$

  The outputed standard deviation of the mean uses the largest possible block averages.

  To improve numerical stability, a slight modification is made to the ensemble average.
  For small volume changes, the ensemble average is very close to unity.
  Thus, the implemented ensemble average is subtracted from unity, but then
  later added back when computing the pressure.
 */
class PressureFromTestVolume : public Modify {
 public:
  //@{
  /** @name Arguments
    - delta_volume: test volume change (default: 1e-4).
   */
  explicit PressureFromTestVolume(argtype args = argtype());
  explicit PressureFromTestVolume(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const override;

  void initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;

  void update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) override;

  std::string write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("PressureFromTestVolume"); }
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<PressureFromTestVolume>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<PressureFromTestVolume>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit PressureFromTestVolume(std::istream& istr);
  ~PressureFromTestVolume();

  //@}
 private:
  double delta_volume_;
};

}  // namespace feasst

#endif  // FEASST_STEPPERS_PRESSURE_FROM_TEST_VOLUME_H_

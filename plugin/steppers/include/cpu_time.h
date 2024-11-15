
#ifndef FEASST_STEPPERS_CPU_TIME_H_
#define FEASST_STEPPERS_CPU_TIME_H_

#include "utils/include/timer.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Periodically print the cpu time spent on the simulation as defined in Timer.
 */
class CPUTime : public AnalyzeWriteOnly {
 public:
  //@{
  /** @name Arguments
    - Stepper arguments.
   */
  explicit CPUTime(argtype args = argtype());
  explicit CPUTime(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void initialize(MonteCarlo * mc) override;
  std::string write(const MonteCarlo& mc) override;

  // serialize
  std::string class_name() const override { return std::string("CPUTime"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<CPUTime>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<CPUTime>(args); }
  explicit CPUTime(std::istream& istr);

  //@}
 private:
  int num_writes_;
  double initialize_time_;
};

inline std::shared_ptr<CPUTime> MakeCPUTime(argtype args = argtype()) {
  return std::make_shared<CPUTime>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CPU_TIME_H_


#ifndef FEASST_STEPPERS_CPU_TIME_H_
#define FEASST_STEPPERS_CPU_TIME_H_

#include "utils/include/timer.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Periodically print the cpu time spent on the simulation.
 */
class CPUTime : public AnalyzeWriteOnly {
 public:
  explicit CPUTime(const argtype& args = argtype());

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  std::string write(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  const Accumulator& accumulator() const override { return steps_per_second_; }

  // serialize
  std::string class_name() const override { return std::string("CPUTime"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<CPUTime>(istr); }
  CPUTime(std::istream& istr);

 private:
  Accumulator steps_per_second_;
  int num_writes_;
  double initialize_time_;
};

inline std::shared_ptr<CPUTime> MakeCPUTime(const argtype &args = argtype()) {
  return std::make_shared<CPUTime>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CPU_TIME_H_

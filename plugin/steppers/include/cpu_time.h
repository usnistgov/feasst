
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
  CPUTime(const argtype& args = argtype());

  void initialize(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    num_writes_ = 0;
    initialize_time_ = cpu_hours();
  }

  std::string write(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    std::stringstream ss;
    ++num_writes_;
    const double elapsed_hours = cpu_hours() - initialize_time_;
    const double steps_per_second = num_writes_*steps_per_write()
                                    /(elapsed_hours*60*60);
    steps_per_second_.accumulate(steps_per_second);
    ss << "steps per second: " << steps_per_second
       << std::endl;
    return ss.str();
  }

  const Accumulator& accumulator() const override { return steps_per_second_; }

  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(1, ostr);
    feasst_serialize(num_writes_, ostr);
    feasst_serialize(initialize_time_, ostr);
  }

  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<CPUTime>(istr); }

  CPUTime(std::istream& istr) : AnalyzeWriteOnly(istr) {
    feasst_deserialize_version(istr);
    feasst_deserialize(&num_writes_, istr);
    feasst_deserialize(&initialize_time_, istr);
  }

  std::string class_name() const override { return std::string("CPUTime"); }

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

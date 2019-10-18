
#ifndef FEASST_STEPPERS_LOG_H_
#define FEASST_STEPPERS_LOG_H_

#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Periodically print the status of the Criteria and Trials, typically in a
  one line format with header.
 */
class Log : public AnalyzeWriteOnly {
 public:
  Log(const argtype& args = argtype());
  void initialize(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    std::stringstream ss;
    ss << criteria->status_header() << " " << trial_factory.status_header()
       << std::endl;
    printer(ss.str());
  }

  std::string write(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    // ensure the following order matches the header from initialization.
    std::stringstream ss;
    ss << criteria->status() << " " << trial_factory.status() << std::endl;
    return ss.str();
  }

  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(1, ostr);
  }

  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Log>(istr); }

  Log(std::istream& istr) : AnalyzeWriteOnly(istr) { feasst_deserialize_version(istr); }

  std::string class_name() const override { return std::string("Log"); }
};

inline std::shared_ptr<Log> MakeLog(const argtype &args = argtype()) {
  return std::make_shared<Log>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_LOG_H_

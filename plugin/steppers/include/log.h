
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
  explicit Log(const argtype& args = argtype());

  void initialize(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  std::string write(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("Log"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Log>(istr); }
  Log(std::istream& istr);
};

inline std::shared_ptr<Log> MakeLog(const argtype &args = argtype()) {
  return std::make_shared<Log>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_LOG_H_


#ifndef FEASST_STEPPERS_LOG_H_
#define FEASST_STEPPERS_LOG_H_

#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Periodically print the status of the Criteria and Trials in a one-line comma
  separated value format with headers.

  This log file is designed more for checking status than computing ensemble
  averages.
  To compute ensemble averages, consider using a custom Analyze or accumulating
  the averages every trial step inside your script.

  The ordering of the columns may change at any time.
  Thus, do not assume a particular order when analyzing the log file.
  Instead, use the headers to specify columns.
  For example, the pandas module in python is ideal for this task.

  By default, the first number printed for a Trial is its Trial::acceptance.
 */
class Log : public AnalyzeWriteOnly {
 public:
  /**
    args:
    - max_precision: use maximum precision if true (default: false).
   */
  explicit Log(argtype args = argtype());
  explicit Log(argtype * args);

  std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const override;

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("Log"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Log>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<Log>(args); }
  Log(std::istream& istr);

 private:
  bool max_precision_;
};

inline std::shared_ptr<Log> MakeLog(argtype args = argtype()) {
  return std::make_shared<Log>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_LOG_H_

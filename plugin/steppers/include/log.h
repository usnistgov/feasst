
#ifndef FEASST_STEPPERS_LOG_H_
#define FEASST_STEPPERS_LOG_H_

#include "system/include/bond_visitor.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Periodically print the status of the Syste, Criteria and Trials in a one-line
  comma separated value format with headers.

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
  //@{
  /** @name Arguments
    - max_precision: use maximum precision if true (default: false).
      See <a href="../../monte_carlo/tutorial/tutorial_0_ref_configs.html">this tutorial</a> for an example.
    - include_bonds: if true, print bond energies (default: true).
    - Stepper arguments.
   */
  explicit Log(argtype args = argtype());
  explicit Log(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  std::string write(const MonteCarlo& mc) override;

  // serialize
  std::string class_name() const override { return std::string("Log"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Log>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<Log>(args); }
  explicit Log(std::istream& istr);

  //@}
 private:
  bool max_precision_;
  bool include_bonds_;
  BondVisitor bond_visitor_;
};

inline std::shared_ptr<Log> MakeLog(argtype args = argtype()) {
  return std::make_shared<Log>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_LOG_H_

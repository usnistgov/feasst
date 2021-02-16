
#ifndef FEASST_UTILS_PROGRESS_REPORT_H_
#define FEASST_UTILS_PROGRESS_REPORT_H_

#include <string>
#include <sstream>
#include <memory>
#include "utils/include/arguments.h"

namespace feasst {

/**
  Periodically report on the progress of long tasks and estimate remaining
  time to completion.
 */
class ProgressReport {
 public:
  /**
    args:
    - num: number of iterations (default: 1).
    - percent_per_write: report every this many percent progress made.
      (default: 0.1 e.g., 10%)
    - file_name: append report to this file. If empty, cout (default: empty).
   */
  explicit ProgressReport(argtype args = argtype());

  /// Set the number of iterations to completion.
  void set_num(const int num) { num_ = num; }

  /// Reset the progress report.
  void reset();

  /// Return the current progress percentage.
  double percent() const;

  /// Check on progress every iteration.
  void check();

  /// Write the report.
  void write();

  /// Serialize object.
  void serialize(std::ostream& ostr) const;

  /// Deserialize object.
  explicit ProgressReport(std::istream& istr);

 private:
  int num_;
  int current_;
  double percent_per_write_;
  double last_percent_ = 0.;
  std::string file_name_;
  double starting_hours_;
};

inline std::shared_ptr<ProgressReport> MakeProgressReport(
    argtype args = argtype()) {
  return std::make_shared<ProgressReport>(args);
}

}  // namespace feasst

#endif  // FEASST_UTILS_PROGRESS_REPORT_H_

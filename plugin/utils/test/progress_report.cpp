#include <fstream>
#include "utils/test/utils.h"
#include "utils/include/progress_report.h"

namespace feasst {

TEST(ProgressReport, percent) {
  int num = 1e6;
  ProgressReport report({{"num", str(num)}, {"percent_per_write", "0.01"},
                         {"file_name", "tmp/prog.txt"},
                         {"double_percent_per_write", "0"}});
  ProgressReport report2({{"num", str(num)}, {"percent_per_write", "1e-4"},
                          {"file_name", "tmp/prog2.txt"},
                          {"double_percent_per_write", "1"}});
  for (int i = 0; i < num; ++i) {
    report.check();
    report2.check();
  }
}

}  // namespace feasst

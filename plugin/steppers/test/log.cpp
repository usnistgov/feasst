#include <fstream>
#include "utils/test/utils.h"
#include "steppers/include/log.h"

namespace feasst {

TEST(Log, serialize) {
  std::ofstream file;
  file.open("templog", std::ofstream::out);
  file << "testing";
  file.close();

  auto log = MakeLog({{"file_name", "templog"},
                      {"clear_file", "true"}});
  auto log2 = test_serialize<Log, Analyze>(*log);
}

}  // namespace feasst

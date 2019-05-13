#include "utils/test/utils.h"
#include "steppers/include/log.h"

namespace feasst {

TEST(Log, serialize) {
  Log log;
  auto log2 = test_serialize<Log, Analyze>(log);
}

}  // namespace feasst

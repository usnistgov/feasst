#include "utils/test/utils.h"
#include "steppers/include/criteria_writer.h"

namespace feasst {

TEST(CriteriaWriter, serialize) {
  CriteriaWriter analyze;
  auto analyze2 = test_serialize<CriteriaWriter, Analyze>(analyze);
}

}  // namespace feasst

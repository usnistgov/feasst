#include "utils/test/utils.h"
#include "flat_histogram/include/write_file_and_check.h"

namespace feasst {

TEST(WriteFileAndCheck, serialize) {
  auto obj = std::make_shared<WriteFileAndCheck>(argtype({{"sim_start", "0"}, {"sim_end", "1"}, {"sim", "0"},
    {"file_prefix", "a"}, {"file_suffix", "a"}, {"output_file", "a"}}));
  Action obj2 = test_serialize(*obj);
}

}  // namespace feasst

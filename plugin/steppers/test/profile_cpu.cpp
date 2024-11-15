#include "utils/test/utils.h"
#include "steppers/include/profile_cpu.h"

namespace feasst {

TEST(ProfileCPU, serialize) {
  auto obj = MakeProfileCPU({{"output_file", "tmp"}});
  auto obj2 = test_serialize_unique(*obj);
}

}  // namespace feasst

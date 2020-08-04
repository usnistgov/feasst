#include "utils/test/utils.h"
#include "utils/include/file.h"

namespace feasst {

TEST(File, exists_and_backup) {
  const std::string name("tmp/tmp345");
  const std::string backup("tmp/tmp345.bak");
  EXPECT_FALSE(file_exists(name));
  std::ofstream file(name);
  file << "hi";
  EXPECT_TRUE(file_exists(name));
  file_backup(name);
  EXPECT_FALSE(file_exists(name));
  EXPECT_TRUE(file_exists(backup));
}

}  // namespace feasst

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

TEST(File, number_lines) {
  EXPECT_EQ(num_lines("../plugin/aniso/test/data/fc.txt"), 6721);
}

TEST(File, compare_files) {
  EXPECT_TRUE (compare_files("../plugin/aniso/test/data/fc.txt", "../plugin/aniso/test/data/fc.txt"));
  EXPECT_FALSE(compare_files("../plugin/aniso/test/data/fc.txt", "../plugin/aniso/test/data/rigid_and_connector.txt"));
}

TEST(File, remove_file) {
  const std::string filename = "tmp/tmp56375";
  std::ofstream file(filename);
  EXPECT_TRUE(file_exists(filename));
  remove_file(filename);
  EXPECT_FALSE(file_exists(filename));
}

}  // namespace feasst

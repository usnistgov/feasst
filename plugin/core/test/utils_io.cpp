#include <gtest/gtest.h>
#include "core/include/utils_io.h"

namespace feasst {

TEST(UtilsIO, split) {
  const std::vector<std::string> str = split("This is a string");
  EXPECT_EQ(4, str.size());
  EXPECT_EQ("This", str[0]);
  EXPECT_EQ("is", str[1]);
  EXPECT_EQ("a", str[2]);
  EXPECT_EQ("string", str[3]);
}

TEST(UtilsIO, trim) {
  const char* f = "/my/path/to/file";
  std::string fs = trim("/", f);
  EXPECT_EQ(0, fs.compare("file"));

  const char* f2 = "file";
  fs.assign(trim("/", f2));
  EXPECT_EQ(0, fs.compare("file"));

  const char* f3 = "/home/username/feasst/forcefield/cg7mabaniso.json";
  EXPECT_EQ("json", trim(".", f3));

  // now find the path by trimming from right (3rd flag 0) instead of the left
  EXPECT_EQ("/home/username/feasst/forcefield/", trim("/", f3, 0));
}

}  // namespace feasst

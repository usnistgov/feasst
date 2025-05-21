#include "utils/test/utils.h"
#include "utils/include/io.h"

namespace feasst {

TEST(IO, split) {
  const std::vector<std::string> str = split("This is a string");
  EXPECT_EQ(4, str.size());
  EXPECT_EQ("This", str[0]);
  EXPECT_EQ("is", str[1]);
  EXPECT_EQ("a", str[2]);
  EXPECT_EQ("string", str[3]);
}

TEST(IO, trim) {
  const char* f = "/my/path/to/file";
  std::string fs = trim("/", f);
  EXPECT_EQ(0, fs.compare("file"));

  const char* f2 = "file";
  fs.assign(trim("/", f2));
  EXPECT_EQ(0, fs.compare("file"));

  const char* f3 = "/home/username/feasst/particle/cg7mabaniso.json";
  EXPECT_EQ("json", trim(".", f3));

  // now find the path by trimming from right (3rd flag 0) instead of the left
  EXPECT_EQ("/home/username/feasst/particle/", trim("/", f3, 0));
}

TEST(IO, serialize) {
  std::vector<double> vec2, vec = {0, 1234.203467642396097823666};
  std::stringstream ss;
  feasst_serialize(vec, ss);
  feasst_deserialize(&vec2, ss);
  EXPECT_EQ(vec, vec2);
}

TEST(IO, serialize_string) {
  std::string str("12bananas");
  std::stringstream ss;
  feasst_serialize(str, ss);
  feasst_serialize(str, ss);
  std::string str2;
  feasst_deserialize(&str2, ss);
  feasst_deserialize(&str2, ss);
  EXPECT_EQ(str, str2);
  EXPECT_EQ(str, "12bananas");
  TRY(
    std::string str3("12 bananas");
    feasst_serialize(str3, ss);
    CATCH_PHRASE("no spaces");
  );
}

TEST(IO, str) {
  std::map<std::string, std::string> mp = {{"hi","you"},{"well","hello"},};
  EXPECT_EQ(str(mp), "hi you well hello ");
}

TEST(IO, is_found_in) {
  EXPECT_TRUE(is_found_in("lets see if you can find this", "this"));
  EXPECT_FALSE(is_found_in("empty", "this"));
}

TEST(IO, sized_int_to_str) {
  EXPECT_EQ("01", sized_int_to_str(1, 16));
}

}  // namespace feasst

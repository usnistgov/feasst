#include "utils/test/utils.h"
#include "utils/include/argument_parse.h"

namespace feasst {

TEST(ArgumentParse, args) {
  char *argv[] = {NULL, (char*)"-k", (char*)"hi", (char*)"-z", (char*)"1234"};
  ArgumentParse args(5, argv);
  EXPECT_EQ(args.str(), "-k,hi,-z,1234,");
  EXPECT_TRUE(args.option_given("-k"));
  EXPECT_FALSE(args.option_given("-f"));
  EXPECT_EQ(args.get("-k"), "hi");
  EXPECT_EQ(args.get("-x", "dflt"), "dflt");
  EXPECT_EQ(args.get_double("-x", 123.5), 123.5);
  EXPECT_EQ(args.get_int("-z"), 1234);
  EXPECT_EQ(args.get_int("-x", -23), -23);
}

}

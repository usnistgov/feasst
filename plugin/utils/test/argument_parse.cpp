#include "utils/test/utils.h"
#include "utils/include/argument_parse.h"

namespace feasst {

TEST(ArgumentParse, args) {
  ArgumentParse args("A test of ArgumentParse", {
    {"-k", "the argument k has a default value", "yo"},
    {"-z", "the argument z has no default"},
    {"-g", "the argument g", "dflt"},
    {"-y", "the argument y", "123.5"},
    {"-q", "the argument q", "-23"},
    {"-p", "the argument p has no default"},
  });
  char *argv[] = {NULL, (char*)"-k", (char*)"hi", (char*)"-z", (char*)"1234"};
  args.parse(5, argv);
  //char *argv[] = {NULL, (char*)"-k", (char*)"hi", (char*)"-z", (char*)"1234",(char*)"-h"};
  //INFO(args.parse(6, argv));
  EXPECT_EQ(args.str(), "-k,hi,-z,1234,-g,dflt,-y,123.5,-q,-23");
  EXPECT_TRUE(args.option_given("-k"));
  EXPECT_FALSE(args.option_given("-f"));
  EXPECT_FALSE(args.option_given("-p"));
  EXPECT_EQ(args.get("-k"), "hi");
  EXPECT_EQ(args.get("-g"), "dflt");
  EXPECT_EQ(args.get_double("-y"), 123.5);
  EXPECT_EQ(args.get_int("-z"), 1234);
  EXPECT_EQ(args.get_int("-q"), -23);
  TRY(
    args.get("--notgiven");
    CATCH_PHRASE("not given");
  );
}

}

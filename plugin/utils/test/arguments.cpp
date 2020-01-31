#include <gtest/gtest.h>
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"

namespace feasst {

TEST(Arguments, Arguments) {
  Arguments args;
  args.init({{"key1", "val1"}});
  EXPECT_EQ(args.key("yo").dflt("hi").str(), "hi");

  // test if a key is empty before it is set
  TRY(
    args.empty();
    CATCH_PHRASE("key must be set before");
  );

  // test if a key is empty after set
  EXPECT_TRUE(args.key("yo").empty());

  // test if all keywords have been used
  TRY(
    EXPECT_FALSE(args.check_all_used());
    CATCH_PHRASE("All keywords provided in args must be used");
  );

  // test if provided key is empty
  EXPECT_FALSE(args.key("key1").empty());

  // check second key pair without providing key, because already set
  EXPECT_EQ(args.dflt("hi").str(), "val1");

  // alternatively, set key and check simultaneously
  EXPECT_EQ(args.key("key1").dflt("hi").str(), "val1");

  // test if key1 is removed using remove()
  const int nargs = args.size();
  EXPECT_EQ(args.key("key1").dflt("hi").remove().str(), "val1");
  EXPECT_EQ(args.size(), nargs - 1);

  // check if all args were used
  EXPECT_TRUE(args.check_all_used());
}

TEST(Arguments, integer) {
  TRY(
    Arguments args;
    args.init({{"key1", "val1"}});
    args.key("key1").integer();
    CATCH_PHRASE("was expected to be an integer");
  );
  TRY(
    Arguments args;
    args.init({{"key1", "1.1"}});
    args.key("key1").integer();
    CATCH_PHRASE("was expected to be an integer");
  );

  Arguments args;
  args.init({{"key1", "1."}});
  args.key("key1").integer();
}

TEST(Arguments, dble) {
  TRY(
    Arguments args;
    args.init({{"key1", "mymypie"}});
    args.key("key1").dble();
    CATCH_PHRASE("was expected to be a double precision floating point number");
  );
  Arguments args;
  args.init({{"key1", "3.1415"}});
  EXPECT_NEAR(3.1415, args.key("key1").dble(), NEAR_ZERO);
}

TEST(Arguments, arglist) {
  arglist argls = {{ {"set1", {{"key1", "val1"}} },
                             {"set2", {{"key2", "val2"}} } }};
  auto set = argls.find("set2");
  auto pair = set->second.find("key2");
  EXPECT_NE(set->second.end(), pair);
  EXPECT_TRUE(pair->second == "val2");
}

TEST(Arguments, boolean) {
  Arguments args(argtype({{"bananas", "1"}}));
  EXPECT_TRUE(args.key("bananas").boolean());
  EXPECT_EQ("used(bananas, )", args.status());
}

}  // namespace feasst

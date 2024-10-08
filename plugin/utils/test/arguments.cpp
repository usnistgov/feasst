#include "utils/test/utils.h"
#include "utils/include/utils.h"
#include "utils/include/arguments_extra.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"

namespace feasst {

class TestArgs {
 public:
  TestArgs(argtype args = argtype()) : TestArgs(&args) {
    feasst_check_all_used(args);
  }
  TestArgs(argtype * args) {
    key1_ = str("strkey", args);
    dblekey_ = dble("dblekey", args);
    intkey_ = integer("intkey", args);
    boolkey_ = boolean("boolkey", args, false);
  }
 private:
  std::string key1_;
  double dblekey_;
  int intkey_;
  bool boolkey_;
};

TEST(Arguments, args) {
  TestArgs({{"strkey", "val1"}, {"dblekey", "0.2"}, {"intkey", "-15"},
            {"boolkey", "true"}});
  TestArgs({{"strkey", "val1"}, {"dblekey", "0.2"}, {"intkey", "-15"}});
  TRY(
    TestArgs({{"strkey", "val1"}, {"dblekey", "0.2"}});
    CATCH_PHRASE("key(intkey) is required for args but not found");
  );
  argtype args;
  EXPECT_EQ(str("yo", &args, "hi"), "hi");
  EXPECT_FALSE(used("yo", args));
  args = {{"key1", "val1"}};
  EXPECT_TRUE(used("key1", args));
  EXPECT_EQ(1, args.size());
  str("key1", &args);
  feasst_check_all_used(args);
}


TEST(Arguments, integer) {
  argtype args = {{"key1", "val1"}};
  TRY(
    integer("key1", &args);
    CATCH_PHRASE("was expected to be an integer");
  );
  args = {{"key1", "1.1"}};
  TRY(
    integer("key1", &args);
    CATCH_PHRASE("was expected to be an integer");
  );
  args = {{"key1", "1"}};
  integer("key1", &args);
}

TEST(Arguments, dble) {
  argtype args = {{"key1", "mymypie"}};
  TRY(
    dble("key1", &args);
    CATCH_PHRASE("was expected to be a double precision number");
  );
  args = {{"key1", "3.1415"}};
  EXPECT_NEAR(3.1415, dble("key1", &args), NEAR_ZERO);
}

//TEST(Arguments, arglist) {
//  arglist argls = {{ {"set1", {{"key1", "val1"}} },
//                             {"set2", {{"key2", "val2"}} } }};
//  auto set = argls.find("set2");
//  auto pair = set->second.find("key2");
//  EXPECT_NE(set->second.end(), pair);
//  EXPECT_TRUE(pair->second == "val2");
//}

TEST(Arguments, boolean) {
  for (const std::string arg : {"1", "true", "True"}) {
    argtype args = {{"bananas", arg}};
    EXPECT_TRUE(boolean("bananas", &args));
  }
}

//TEST(Arguments, parse) {
//  arglist args = {{{"Random", {{"seed", "time"}}}}};
//  EXPECT_EQ(1, args.size());
//  argtype random_args = get("Random", &args);
//  EXPECT_EQ(0, args.size());
//  EXPECT_EQ(1, random_args.size());
//  EXPECT_EQ("time", str("seed", &random_args));
//  EXPECT_EQ(0, random_args.size());
//}

TEST(Arguments, str) {
  argtype arg = {{"hi", "you"}};
  EXPECT_EQ("hi you ", str(arg));
  arglist args = {{{"major_key1", {{"minor_key1", "value1"}}}}};
  EXPECT_EQ("{{{\"major_key1\",minor_key1 value1 },}}", str(args));
}

TEST(Utils, find_in_list) {
  arglist args = {{
    {"major_key", {{"minor_key", "value1"}}},
    {"major_key2", {{"minor_key2", "value2"}, {"minor_key3", "value3"}}}
  }};
  int find;
  find_in_list(std::string("major_key2"), args, &find);
  EXPECT_EQ(1, find);
  // INFO(str(args));
  replace_value("value3", "value[sim_index]", &args);
  EXPECT_EQ(str(args), "{{{\"major_key\",minor_key value1 },{\"major_key2\",minor_key2 value2 minor_key3 value[sim_index] },}}");
  replace_in_value("[sim_index]", "5", &args);
  EXPECT_EQ(str(args), "{{{\"major_key\",minor_key value1 },{\"major_key2\",minor_key2 value2 minor_key3 value5 },}}");
}

TEST(Arguments, line_to_argtype) {
  argtype args = line_to_argtype("key1 val1 key2 val2");
  argtype expected = {{"key1", "val1"}, {"key2", "val2"}};
  EXPECT_EQ(args, expected);
}

TEST(Arguments, parse_dimensional) {
  argtype args = {{"x0", "1"}, {"x1", "2"}};
  std::vector<double> dat = parse_dimensional("x", &args, 4);
  EXPECT_EQ(args.size(), 0);
  EXPECT_EQ(dat.size(), 2);
  EXPECT_EQ(dat[0], 1);
  EXPECT_EQ(dat[1], 2);
}

}  // namespace feasst

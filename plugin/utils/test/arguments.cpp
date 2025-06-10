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
  EXPECT_EQ("hi=you", str(arg));
  arglist args = {{{"major_key1", {{"minor_key1", "value1"}}}}};
  EXPECT_EQ("{{{\"major_key1\",minor_key1=value1},}}", str(args));
}

TEST(Utils, find_in_list) {
  arglist args = {{
    {"major_key", {{"minor_key", "value1"}}},
    {"major_key2", {{"minor_key2", "value2"}, {"minor_key3", "value3"}}}
  }};
  int find;
  find_in_list(std::string("major_key2"), args, &find);
  EXPECT_EQ(1, find);
  // DEBUG(str(args));
  replace_value("value3", "value[sim_index]", &args);
  EXPECT_EQ(str(args), "{{{\"major_key\",minor_key=value1},{\"major_key2\",minor_key2=value2 minor_key3=value[sim_index]},}}");
  replace_in_value("[sim_index]", "5", &args);
  EXPECT_EQ(str(args), "{{{\"major_key\",minor_key=value1},{\"major_key2\",minor_key2=value2 minor_key3=value5},}}");
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

TEST(Arguments, parse_line) {
  std::pair<std::string, argtype> parsed;
  // deprecated old interface style
  parsed = parse_line("Configuration cubic_box_length 12 particle_type lj.txt ", NULL, NULL);
  EXPECT_EQ("Configuration", parsed.first);
  EXPECT_EQ(parsed.second.find("cubic_box_length")->first, "cubic_box_length");
  EXPECT_EQ(parsed.second.find("cubic_box_length")->second, "12");
  EXPECT_EQ(parsed.second.find("particle_type")->first, "particle_type");
  EXPECT_EQ(parsed.second.find("particle_type")->second, "lj.txt");

  parsed = parse_line("Configuration cubic_box_length=12=5 particle_type=lj.txt ", NULL, NULL);
  EXPECT_EQ("Configuration", parsed.first);
  EXPECT_EQ(parsed.second.find("cubic_box_length")->first, "cubic_box_length");
  EXPECT_EQ(parsed.second.find("cubic_box_length")->second, "12=5");
  EXPECT_EQ(parsed.second.find("particle_type")->first, "particle_type");
  EXPECT_EQ(parsed.second.find("particle_type")->second, "lj.txt");
}

TEST(Arguments, parse_mcs_for) {
  std::stringstream ss;
  ss << "For [len]=8,9\nFor [pt]=lj,atom\nConfiguration cubic_box_length [len] particle_type [pt].txt\nEndFor\nEndFor";
  std::vector<arglist> mcs = parse_mcs(ss);
  for (const auto& lp : mcs[0]) {
    EXPECT_EQ(lp.first, "Configuration");
  }
  EXPECT_EQ(mcs[0][0].second["cubic_box_length"], "8");
  EXPECT_EQ(mcs[0][1].second["cubic_box_length"], "8");
  EXPECT_EQ(mcs[0][2].second["cubic_box_length"], "9");
  EXPECT_EQ(mcs[0][3].second["cubic_box_length"], "9");
  EXPECT_EQ(mcs[0][0].second["particle_type"], "lj.txt");
  EXPECT_EQ(mcs[0][1].second["particle_type"], "atom.txt");
  EXPECT_EQ(mcs[0][2].second["particle_type"], "lj.txt");
  EXPECT_EQ(mcs[0][3].second["particle_type"], "atom.txt");
  //DEBUG(mcs[0][0].first << " : " << feasst_str(mcs[0][0].second));
  TRY(
    std::stringstream ss2;
    ss2 << "For [obj]=Configuration,8";
    parse_mcs(ss2);
    CATCH_PHRASE("Found a For without a corresponding EndFor");
  );
  TRY(
    std::stringstream ss2;
    ss2 << "For [var1]:[var2]=0:1,2:3,4\nEndFor";
    parse_mcs(ss2);
    CATCH_PHRASE("did not equal the number of given");
  );
  TRY(
    std::stringstream ss2;
    ss2 << "For var1:[var2]=0:1,2:4\nEndFor";
    parse_mcs(ss2);
    CATCH_PHRASE("should begin with a ");
  );
  TRY(
    std::stringstream ss2;
    ss2 << "For [var1:[var2]=0:1,2:4\nEndFor";
    parse_mcs(ss2);
    CATCH_PHRASE("should end with a ");
  );
}

TEST(Arguments, parse_mcs_for_if) {
  std::stringstream ss;
  ss << "For [len]=8,9\nIf defined=?hi\nRemove\nElse\nConfiguration cubic_box_length [len] particle_type lj.txt\nEndIf\nEndFor";
  std::vector<arglist> mcs = parse_mcs(ss);
  for (const auto& lp : mcs[0]) {
    EXPECT_EQ(lp.first, "Remove");
  }
  std::stringstream ss2;
  ss2 << "For [len]=8,9\nIf defined=?\nRemove\nElse\nConfiguration cubic_box_length [len] particle_type lj.txt\nEndIf\nEndFor";
  mcs = parse_mcs(ss2);
  for (const auto& lp : mcs[0]) {
    EXPECT_EQ(lp.first, "Configuration");
  }
  std::stringstream ss3;
  ss3 << "For [len]:[pt]=8:?8,9:?\nConf [len]=[pt]\nEndFor";
  mcs = parse_mcs(ss3);
  EXPECT_EQ(static_cast<int>(mcs.size()), 1);
  EXPECT_EQ(static_cast<int>(mcs[0].size()), 2);
  EXPECT_EQ(mcs[0][0].first, "Conf");
  EXPECT_EQ(str(mcs[0][0].second), "8=8");
  EXPECT_EQ(mcs[0][1].first, "Conf");
  EXPECT_EQ(str(mcs[0][1].second), "9=?");
  std::stringstream ss4;
  ss4 << "For [cutoff]:[xyz]=28.741260651153034:?,14:?\n\
Let [Config]=Configuration particle_type=trappe:/feasst/particle/ethane.txt cutoff=[cutoff]\n\
If defined=?[xyz]\n\
[Config] xyz_file=?[xyz]\n\
Else\n\
[Config] cubic_side_length=?71.85315162788258\n\
Endif\n\
Endfor";
  mcs = parse_mcs(ss4);
}

TEST(Arguments, IF) {
  TRY(
    std::stringstream ss;
    ss << "If something=something";
    parse_mcs(ss);
    CATCH_PHRASE("syntax is");
  );
  TRY(
    std::stringstream ss;
    ss << "If defined=?";
    parse_mcs(ss);
    CATCH_PHRASE("If without a corresponding EndIf");
  );
  TRY(
    std::stringstream ss;
    ss << "If defined=?\nEndIf\nEndIf";
    parse_mcs(ss);
    CATCH_PHRASE("EndIf without a corresponding If");
  );
  TRY(
    std::stringstream ss;
    ss << "Else\nEndIf";
    parse_mcs(ss);
    CATCH_PHRASE("Else without a corresponding If");
  );
  std::stringstream ss;
  ss << "If defined=?\nNot\nElse\nTrue\nEndIf";
  std::vector<arglist> list = parse_mcs(ss);
  EXPECT_EQ("True", list[0][0].first);
  std::stringstream ss2;
  ss2 << "If defined=?Yes\nNot\nElse\nTrue\nEndIf";
  list = parse_mcs(ss2);
  EXPECT_EQ("Not", list[0][0].first);
  std::stringstream ss3;
  ss3 << "If undefined=?Yes\nNot\nElse\nTrue\nEndIf";
  list = parse_mcs(ss3);
  EXPECT_EQ("True", list[0][0].first);
}

}  // namespace feasst

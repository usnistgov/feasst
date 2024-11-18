#include <limits>
#include <sstream>
#include "utils/test/utils.h"
#include "utils/include/serialize.h"

namespace feasst {

typedef std::vector<std::vector<std::vector<double> > > vec3;

TEST(Serialize, serialize) {
  vec3 data = { { {1, 1}, {2, 2} }, { {3, 3}, {4, 4} } };
  std::stringstream ss;
  feasst_serialize(data, ss);
  vec3 data2;
  feasst_deserialize(&data2, ss);
  EXPECT_EQ(data, data2);
}

TEST(Serialize, inf) {
  const long double inf = 2*std::numeric_limits<long double>::max();
  std::stringstream ss;
  feasst_serialize(inf, ss);
  DEBUG(ss.str());
  long double inf2 = -1;
  feasst_deserialize(&inf2, ss);
  DEBUG(inf2);
  EXPECT_EQ(inf, inf2);
}

TEST(Serialize, zero) {
  const double zero = 1e-320;
  std::stringstream ss;
  feasst_serialize(zero, ss);
  DEBUG(ss.str());
  double zero2;
  feasst_deserialize(&zero2, ss);
  DEBUG(zero2);
  EXPECT_EQ(zero, zero);

  const double z = 0;
  feasst_serialize(z, ss);
  //INFO(ss.str());
}

TEST(Serialize, argtype) {
  const argtype& args = {{"hi", "you"}, {"0", "12"}};
  std::stringstream ss;
  feasst_serialize(args, ss);
  EXPECT_EQ(ss.str(), "2 1 0 1 12 1 hi 1 you ");
  argtype args2;
  feasst_deserialize(&args2, ss);
  EXPECT_EQ(args, args2);
  EXPECT_EQ(args2["hi"], "you");
  EXPECT_EQ(args2["0"], "12");
}

TEST(Serialize, arglist) {
  const arglist& args = {{{"hi", {{"0", "12"}}}}};
  std::stringstream ss;
  feasst_serialize(args, ss);
  EXPECT_EQ(ss.str(), "1 1 hi 1 1 0 1 12 ");
  arglist args2;
  feasst_deserialize(&args2, ss);
  //INFO(str(args2));
  EXPECT_EQ(args2[0].first, "hi");
  EXPECT_EQ(args2[0].second["0"], "12");
}

}  // namespace feasst

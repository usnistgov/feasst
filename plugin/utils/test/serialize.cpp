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

}  // namespace feasst

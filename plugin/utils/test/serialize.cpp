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

}  // namespace feasst

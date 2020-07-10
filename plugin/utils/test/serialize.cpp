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

}  // namespace feasst

#include <sstream>
#include <gtest/gtest.h>
#include "core/include/movie.h"

namespace feasst {

TEST(Movie, serialize) {
  Movie analyze;
  std::stringstream ss, ss2;
  analyze.serialize(ss);
  EXPECT_EQ("Movie 1 1 0 0 1 ", ss.str());
  std::shared_ptr<Analyze> analyze2 = Movie().deserialize(ss);
  analyze2->serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst

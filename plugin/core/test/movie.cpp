#include <sstream>
#include <gtest/gtest.h>
#include "core/include/movie.h"

namespace feasst {

TEST(Movie, serialize) {
  auto movie = MakeMovie({{"file_name", "tmp"}});
  std::stringstream ss, ss2;
  movie->serialize(ss);
  std::shared_ptr<Analyze> movie2 = movie->deserialize(ss);
  movie2->serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst

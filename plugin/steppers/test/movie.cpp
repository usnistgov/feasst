#include "utils/test/utils.h"
#include "steppers/include/movie.h"

namespace feasst {

TEST(Movie, serialize) {
  auto movie = MakeMovie({{"file_name", "tmp"}});
  auto movie2 = test_serialize<Movie, Analyze>(*movie);
}

}  // namespace feasst

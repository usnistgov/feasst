#include "utils/test/utils.h"
#include "math/include/utils_math.h"
#include "system/include/hard_sphere.h"
#include "monte_carlo/include/monte_carlo.h"
#include "models/include/square_well.h"
#include "mayer/include/mayer_sampling.h"
#include "patch/include/visit_model_inner_patch.h"
#include "patch/include/file_xyz_patch.h"
#include "patch/include/movie_patch.h"

namespace feasst {

TEST(MoviePatch, serialize) {
  auto patch = MakeMoviePatch({{"output_file", "hi"}});
  MoviePatch patch2 = test_serialize(*patch);
}

}  // namespace feasst

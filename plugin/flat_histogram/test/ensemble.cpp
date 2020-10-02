#include "utils/test/utils.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/ensemble.h"

namespace feasst {

TEST(Ensemble, reweight) {
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", "5"}, {"min", "0"}})),
    MakeTransitionMatrix({{"min_sweeps", "10"}}),
    {{"beta", str(1./1.5)},
     {"chemical_potential", "-2.352321"}});
//    {{"soft_max", "5"}, {"soft_min", "1"}}));
  LnProbability lnpirw = Ensemble(*criteria).reweight(1.5);
  EXPECT_NEAR(-7.75236, lnpirw.values()[0], 1e-5);
}

}  // namespace feasst

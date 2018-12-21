#include <gtest/gtest.h>
#include "core/include/visit_model_intra.h"
#include "core/include/model_lj.h"
#include "core/test/configuration_test.h"

namespace feasst {

TEST(VisitModelIntra, energy) {
  Configuration config = chain10_sample();
  ModelLJ model;
  VisitModelIntra visit;
  // don't compute intraparticle interactions between bonded sites.
  visit.set_intra_cut(1);
  // set cut-off to 2.5 so only beads 2 away can interact.
  // due to periodic boundary conditions matching exactly the length,
  // each bead interacts twice at a distance of 2
  config.set_model_param("cutoff", 0, 2.5);
  model.compute(visit, config);
  EXPECT_NEAR(10*2*(4*(pow(2, -12)-pow(2, -6))), visit.energy(), NEAR_ZERO);
}

}  // namespace feasst

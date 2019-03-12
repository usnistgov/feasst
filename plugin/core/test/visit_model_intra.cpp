#include <gtest/gtest.h>
#include "core/include/visit_model_intra.h"
#include "core/include/model_lj.h"
#include "core/test/configuration_test.h"

namespace feasst {

TEST(VisitModelIntra, energy) {
  Configuration config = chain10_sample();
  ModelLJ model;
  VisitModelIntra visit;
  // set cut-off to 2.5 so only beads 2 away can interact.
  // due to periodic boundary conditions matching exactly the length,
  // each bead interacts once each at a distance of 2
  config.set_model_param("cutoff", 0, 2.5);

  // don't compute intraparticle interactions between bonded sites.
  visit.set_intra_cut(1);

  model.compute(&config, &visit);
  const double pe_lj =  4*(pow(2, -12)-pow(2, -6));
  EXPECT_NEAR(10*pe_lj, visit.energy(), NEAR_ZERO);

  Select all = config.selection_of_all();
  Select ignore;
  ignore.add_site(0, 0);
  all.exclude(ignore);

  model.compute(all, &config, &visit);
  EXPECT_NEAR(8*pe_lj, visit.energy(), NEAR_ZERO);
}

}  // namespace feasst

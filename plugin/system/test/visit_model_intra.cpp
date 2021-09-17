#include <cmath>
#include "utils/test/utils.h"
#include "configuration/include/domain.h"
#include "configuration/include/utils.h"
#include "system/include/visit_model_intra.h"
#include "system/include/lennard_jones.h"

namespace feasst {

TEST(VisitModelIntra, energy) {
  Configuration config(MakeDomain({{"cubic_box_length", "10"}}),
    {{"particle_type", "../forcefield/chain10.fstprt"}});
  config.add_particle_of_type(0);
  LennardJones model;
  // don't compute intraparticle interactions between bonded sites.
  auto visit = MakeVisitModelIntra({{"cutoff", "1"}});
  visit->precompute(&config);
  // set cut-off to 2.5 so only beads 2 away can interact.
  // due to periodic boundary conditions matching exactly the length,
  // each bead interacts once each at a distance of 2
  config.set_model_param("cutoff", 0, 2.5);

  model.compute(&config, visit.get());
  const double pe_lj =  4*(std::pow(2, -12) - std::pow(2, -6));
  EXPECT_NEAR(10*pe_lj, visit->energy(), NEAR_ZERO);

  // test exclusion of first site (no 0-2 or 0-8 interaction)
  Select all = config.selection_of_all();
  Select ignore;
  ignore.add_site(0, 0);
  all.exclude(ignore);

  model.compute(all, &config, visit.get());
  EXPECT_NEAR(8*pe_lj, visit->energy(), NEAR_ZERO);

  auto visit2 = test_serialize<VisitModelIntra, VisitModel>(*visit);
  model.compute(all, &config, visit2.get());
  EXPECT_NEAR(visit->energy(), visit2->energy(), NEAR_ZERO);
}

}  // namespace feasst

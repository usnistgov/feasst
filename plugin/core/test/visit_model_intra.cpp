#include <gtest/gtest.h>
#include "core/include/visit_model_intra.h"
#include "core/include/model_lj.h"

TEST(VisitModelIntra, energy) {
  feasst::Configuration config;
  config.set_domain(feasst::Domain().set_cubic(10));
  config.add_particle_type("../forcefield/data.chain10");
  // set cut-off to 2.5 so only beads 2 away can interact.
  config.set_model_param("cutoff", 0, 2.5);
  config.add_particle(0);
  feasst::ModelLJ model;
  feasst::VisitModelIntra visit;
  // don't compute intraparticle interactions between bonded sites.
  visit.set_intra_cut(1);
  EXPECT_EQ(config.num_particles(), 1);
  model.compute(visit, config, config.selection_of_all());
  // due to periodic boundary conditions matching exactly the length,
  // each bead interacts twice at a distance of 2
  EXPECT_NEAR(10*2*(4*(pow(2, -12)-pow(2, -6))), visit.energy(), feasst::NEAR_ZERO);
}


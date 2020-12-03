#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "math/include/constants.h"
#include "configuration/include/file_lmp.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/utils.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model.h"
#include "system/include/model_empty.h"
#include "system/include/model_two_body_factory.h"
#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

TEST(VisitModel, spce_reference_config) {
  Configuration config = spce_sample1();
  LennardJones model;
  VisitModel visit;
  visit.precompute(&config);
  visit.compute(&model, &config);
  const double pe_lj = 827.61105444941393;
  EXPECT_NEAR(pe_lj, visit.energy(), feasst::NEAR_ZERO);
  ModelEmpty empty;
  LongRangeCorrections lrc;
  empty.compute(&config, &lrc);
  const double pe_lrc = -6.8487471455514575;
  EXPECT_NEAR(pe_lrc, lrc.energy(), 1e-13);

//  // test adding/deleting particles, resulting in a ghost
//  Select select;
//  TrialSelectParticle tselect;
//  RandomMT19937 random;
//  tselect.random_particle(config, &select, &random);
//  visit.compute(model, select, &config);
//  const double pe_previous = visit.energy();
//  const Particle& spart = config.select_particle(select.particle_index(0));
//  const double x1_previous = spart.site(0).position().coord(1);
//  config.add_particle_of_type(0);
//  Select new_part;
//  new_part.last_particle_added(&config);
//  config.replace_position(new_part, spart);
//  config.remove_particle(select);
//  visit.compute(model, &config);
//  EXPECT_NEAR(pe_lj, visit.energy(), 1e-12);
//  empty.compute(&config, &lrc);
//  EXPECT_NEAR(pe_lrc, lrc.energy(), 1e-13);
//  EXPECT_EQ(101, config.particles().num()); // includes one ghost particle
//  EXPECT_EQ(100, config.selection_of_all().num_particles());
//  config.check();
//  visit.compute(model, new_part, &config);
//  EXPECT_NEAR(pe_previous, visit.energy(), NEAR_ZERO);
//  const double x2_previous = config.select_particle(new_part.particle_index(0)).site(0).position().coord(1);
//  EXPECT_EQ(x1_previous, x2_previous);
//  visit.check_energy(model, &config);

  // check that the energy of the deleted ghost is the same as the new replacement
  // HWH config.remove_particle(new_
}

}  // namespace feasst

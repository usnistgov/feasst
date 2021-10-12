#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "math/include/constants.h"
#include "configuration/include/file_lmp.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/select.h"
#include "configuration/test/config_utils.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model.h"
#include "system/include/model_empty.h"
#include "system/include/model_two_body_factory.h"

namespace feasst {

double en_lj(const double pos) {
  return 4*(pow(pos, -12) - pow(pos, -6));
}

TEST(VisitModel, energy) {
  Configuration config = two_particle_configuration();
  const double pos = 1.25;
  EXPECT_EQ(config.particle(0).site(0).position().coord(0), 0);
  EXPECT_EQ(config.particle(1).site(0).position().coord(0), pos);
  LennardJones model;
  VisitModel visit;
  visit.precompute(&config);
  visit.compute(&model, &config);
  EXPECT_NEAR(en_lj(pos), visit.energy(), NEAR_ZERO);

  // check PBCs
  Position position = config.particle(1).site(0).position();
  position.set_coord(0, 4);

  Particle particle = config.particle(1);
  Site site = particle.site(0);
  site.set_position(position);
  particle.set_site(0, site);
  Select select;
  select.add_particle(config.select_particle(1), 1);
  config.replace_position(select, particle);
  EXPECT_EQ(3, config.particle(0).site(0).position().size());
  EXPECT_EQ(0, config.particle(0).site(0).position().coord(0));
  EXPECT_EQ(4, config.particle(1).site(0).position().coord(0));
  model.compute(select, &config, &visit);
  EXPECT_NEAR(en_lj(2.), visit.energy(), NEAR_ZERO);
  Select select2;
  select2.add_particle(config.select_particle(0), 0);
  model.compute(select2, &config, &visit);
  EXPECT_NEAR(en_lj(2.), visit.energy(), NEAR_ZERO);

  // serialize
  auto visit2 = test_serialize<VisitModel, VisitModel>(visit);
  EXPECT_NEAR(en_lj(2.), visit2->energy(), NEAR_ZERO);
}

TEST(VisitModel, reference_config) {
  Configuration config = lj_sample4();
  LennardJones model;
  VisitModel visit;
  visit.precompute(&config);
  model.compute(&config, &visit);
  const double en_lj_expect = -16.790321304625856;
  EXPECT_NEAR(en_lj_expect, visit.energy(), NEAR_ZERO);
  const double energy_prev = visit.energy();
  ModelEmpty empty;
  LongRangeCorrections lrc;
  empty.compute(&config, &lrc);
  EXPECT_NEAR(-0.5451660014945704, lrc.energy(), NEAR_ZERO);
  visit.check_energy(&model, &config);

  // test factory double counts with two identical LJ models.
  ModelTwoBodyFactory factory({MakeLennardJones(), MakeLennardJones()});
  factory.compute(&config, &visit);
  EXPECT_NEAR(2.*energy_prev, visit.energy(), NEAR_ZERO);

  // Energy map is not used by default
  EXPECT_FALSE(visit.inner().is_energy_map());

  // Find energy of one particle
  Select one;
  one.add_particle(config.particle(0), 0);
  visit.compute(&model, one, &config);
  EXPECT_NEAR(-3.2639025245521616, visit.energy(), NEAR_ZERO);
}

}  // namespace feasst

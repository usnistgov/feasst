#include <gtest/gtest.h>
#include <random>
#include "core/include/perturb_translate.h"
#include "core/include/perturb_transfer.h"
#include "core/include/model_lj.h"
#include "core/test/system_test.h"
#include "core/test/particle_test.h"

namespace feasst {

TEST(Perturb, Revert) {
  seed_random_by_date();
  System system = default_system();
  PerturbTransfer attempt;
  Particle particle = default_particle();
  particle.set_coordinate(1, 1.25);
  Site site = particle.site(0);
  site.set_coordinate(1, 1.25);
  particle.set_site(0, site);
  System system2(system);
  ModelLJ model;
  VisitModel visit;
  Configuration* config = system.get_configuration();
  model.compute(config, &visit);
  const double pe_original = 4*(pow(1.25, -12) - pow(1.25, -6));
  EXPECT_NEAR(pe_original, visit.energy(), NEAR_ZERO);
  attempt.add(particle, &system);
  EXPECT_EQ(config->num_particles(), 3);
  model.compute(config, &visit);
  const double tri_distance = sqrt(1.25*1.25 + 1.25*1.25);
  const double peTri = 4*(pow(tri_distance, -12) - pow(tri_distance, -6));
  EXPECT_NEAR(2*pe_original + peTri, visit.energy(), NEAR_ZERO);
  attempt.revert();
  EXPECT_EQ(config->num_particles(), 2);
  model.compute(config, &visit);
  EXPECT_NEAR(pe_original, visit.energy(), NEAR_ZERO);

  model.compute(config, &visit);
  EXPECT_NEAR(pe_original, visit.energy(), NEAR_ZERO);
  Configuration* config2 = system2.get_configuration();
  //config2->select_particle(1);
  attempt.select_random_particle(0, config);
  attempt.remove_selected_particle(&system2);
  EXPECT_EQ(1, config2->num_particles());
  model.compute(config2, &visit);
  EXPECT_EQ(0., visit.energy());
  attempt.revert();
  EXPECT_EQ(2, config2->num_particles());
  model.compute(config, &visit);
  EXPECT_NEAR(pe_original, visit.energy(), NEAR_ZERO);

  Position trajectory({0., 0., 0.});
  // trajectory.set_to_origin_3D();
  trajectory.set_coord(2, 1.25);
  PerturbTranslate attempt2;
  attempt2.select_random_particle(0, *config);
  attempt2.translate_selected_particle(trajectory, &system2);
}

}  // namespace feasst

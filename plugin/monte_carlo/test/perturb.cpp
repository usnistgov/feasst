#include <gtest/gtest.h>
#include <random>
#include "monte_carlo/include/perturb.h"
#include "system/include/model_lj.h"
#include "system/test/system_test.h"
#include "configuration/test/particle_test.h"

namespace feasst {

TEST(Perturb, Revert) {
  seed_random_by_date();
  System system = default_system();
  PerturbAdd add;
  System system2(system);
  ModelLJ model;
  VisitModel visit;
  Configuration* config = system.get_configuration();
  model.compute(config, &visit);
  const double pe_original = 4*(pow(1.25, -12) - pow(1.25, -6));
  EXPECT_NEAR(pe_original, visit.energy(), NEAR_ZERO);
  Position position;
  position.set_vector({0, 1.25, 0});
  auto tsel_ghost = MakeTrialSelectParticle({{"particle_type", "0"}});
  // precompute adds argument {"ghost", "true"} to tsel
  add.precompute(tsel_ghost.get(), &system);
  tsel_ghost->select(Select(), &system);
  add.add(&system, tsel_ghost.get(), position);
  EXPECT_EQ(config->num_particles(), 3);
  model.compute(config, &visit);
  const double tri_distance = sqrt(1.25*1.25 + 1.25*1.25);
  const double peTri = 4*(pow(tri_distance, -12) - pow(tri_distance, -6));
  EXPECT_NEAR(2*pe_original + peTri, visit.energy(), NEAR_ZERO);
  add.revert(&system);
  EXPECT_EQ(config->num_particles(), 2);
  model.compute(config, &visit);
  EXPECT_NEAR(pe_original, visit.energy(), NEAR_ZERO);

  model.compute(config, &visit);
  EXPECT_NEAR(pe_original, visit.energy(), NEAR_ZERO);
  Configuration* config2 = system2.get_configuration();
  PerturbRemove remove;
  auto tsel = MakeTrialSelectParticle({{"particle_type", "0"}});
  tsel->select(Select(), &system2);
  remove.perturb(&system2, tsel.get());
  EXPECT_EQ(2, config2->num_particles());
  remove.finalize(&system2);
  EXPECT_EQ(1, config2->num_particles());
  model.compute(config2, &visit);
  EXPECT_EQ(0., visit.energy());

  Position trajectory({0., 0., 0.});
  trajectory.set_coord(2, 1.25);
  PerturbTranslate attempt2;
  tsel->select(Select(), &system);
  attempt2.perturb(&system, tsel.get());
}

}  // namespace feasst

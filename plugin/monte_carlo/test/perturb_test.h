#include <gtest/gtest.h>
#include <random>
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "system/include/lennard_jones.h"
#include "system/test/system_test.h"
#include "configuration/test/particle_test.h"

namespace feasst {

// Test the ability to revert addition, removal and translation
// Make this a function so other tests can call it with different systems.
inline void test_revert(System * system) {
  const double pe_original = 4*(pow(1.25, -12) - pow(1.25, -6));
  EXPECT_NEAR(pe_original, system->energy(), NEAR_ZERO);

  const VisitModelInner * inner = system->potential(0).visit_model()->inner();
  system->finalize();
  if (inner->is_energy_map_queryable()) {
    EXPECT_NEAR(pe_original, inner->energy_map()->total_energy(), 10*NEAR_ZERO);
  }
  system->check();

  // translate a particle and revert
//  Position trajectory({0., 0., 0.});
//  trajectory.set_coord(2, 1.25);
  PerturbTranslate translate;
  auto tsel = MakeTrialSelectParticle({{"particle_type", "0"}});
  auto random = MakeRandomMT19937();
//  random = MakeRandomMT19937({{"seed", "default"}});
//  random = MakeRandomMT19937({{"seed", "1580853628"}});
  tsel->select(Select(), system, random.get());
  translate.perturb(system, tsel.get(), random.get());
  EXPECT_GT(std::abs(pe_original - system->energy()), NEAR_ZERO);
//  if (inner->is_energy_map_queryable()) {
//    INFO(inner->energy_map()->total_energy());
//    EXPECT_GT(std::abs(pe_original - inner->energy_map()->total_energy()), NEAR_ZERO);
//  }
  translate.revert(system);
  EXPECT_NEAR(pe_original, system->energy(), NEAR_ZERO);
  if (inner->is_energy_map_queryable()) {
    EXPECT_NEAR(pe_original, inner->energy_map()->total_energy(), 10*NEAR_ZERO);
  }
  system->check();

  // add a third particle, but revert
  Position position;
  position.set_vector({0, 1.25, 0});
  auto tsel_ghost = MakeTrialSelectParticle({{"particle_type", "0"}});
  // precompute adds argument {"ghost", "true"} to tsel
  PerturbAdd add;
  add.precompute(tsel_ghost.get(), system);
  tsel_ghost->select(Select(), system, random.get());
  add.add(system, tsel_ghost.get(), random.get(), position);
  EXPECT_EQ(system->configuration().num_particles(), 3);
  const double tri_distance = sqrt(1.25*1.25 + 1.25*1.25);
  const double peTri = 4*(pow(tri_distance, -12) - pow(tri_distance, -6));
  EXPECT_NEAR(2*pe_original + peTri, system->energy(), NEAR_ZERO);

//  if (inner->is_energy_map_queryable()) {
//    EXPECT_NEAR(2.*pe_original + peTri, inner->energy_map()->total_energy(), 10*NEAR_ZERO);
//  }

  // revert the addition of the third particle
  add.revert(system);
  EXPECT_EQ(system->configuration().num_particles(), 2);
  EXPECT_NEAR(pe_original, system->energy(), NEAR_ZERO);
  if (inner->is_energy_map_queryable()) {
    EXPECT_NEAR(pe_original, inner->energy_map()->total_energy(), 10*NEAR_ZERO);
  }
  system->check();

  // remove one of the remaining two particles
  PerturbRemove remove;
  tsel->select(Select(), system, random.get());
  remove.perturb(system, tsel.get(), random.get());
  EXPECT_EQ(2, system->configuration().num_particles());

  // revert removal
  remove.revert(system);

  EXPECT_NEAR(pe_original, system->energy(), NEAR_ZERO);
  if (inner->is_energy_map_queryable()) {
    EXPECT_NEAR(pe_original, inner->energy_map()->total_energy(), 10*NEAR_ZERO);
  }
  system->check();

  // finalize removal
  remove.perturb(system, tsel.get(), random.get());
  remove.finalize(system);
  EXPECT_EQ(1, system->configuration().num_particles());
  if (inner->is_energy_map_queryable()) {
    EXPECT_NEAR(0., inner->energy_map()->total_energy(), 10*NEAR_ZERO);
  }
  system->check();
  EXPECT_EQ(0., system->energy());

  // add particle (finalize)
  tsel_ghost->select(Select(), system, random.get());
  add.perturb(system, tsel_ghost.get(), random.get());
  const double en2 = system->energy();
  add.finalize(system);
  EXPECT_EQ(2, system->configuration().num_particles());
  EXPECT_GT(std::abs(en2 - pe_original), NEAR_ZERO);
//  INFO("en2 " << en2);
  if (inner->is_energy_map_queryable()) {
    EXPECT_NEAR(en2, inner->energy_map()->total_energy(), 10*NEAR_ZERO);
  }
  system->check();

  // translate (finalize)
  tsel->select(Select(), system, random.get());
  EXPECT_EQ(2, system->configuration().num_particles());
  translate.perturb(system, tsel.get(), random.get());
  const double en2t = system->energy();
  translate.finalize(system);
  if (en2t != 0 || en2 != 0) {
    EXPECT_GT(std::abs(en2t - en2), NEAR_ZERO);
    //INFO("en2t " << en2t);
  }
  system->check();
  //INFO(system->energy());
}

}  // namespace feasst

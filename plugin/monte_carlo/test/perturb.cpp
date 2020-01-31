#include <gtest/gtest.h>
#include <random>
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "system/include/lennard_jones.h"
#include "system/test/system_test.h"
#include "configuration/test/particle_test.h"
#include "system/include/energy_map_all.h"

namespace feasst {

// Test the ability to revert addition, removal and translation
TEST(Perturb, Revert) {
  // try both without and with map
  for (int map = 0; map <= 1; ++map) {
    INFO("map: " << map);
    System system = default_system();
    if (map == 1) {
      system.set_unoptimized(0,
        Potential(MakeLennardJones(),
                  MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
    }
    const double pe_original = 4*(pow(1.25, -12) - pow(1.25, -6));
    EXPECT_NEAR(pe_original, system.energy(), NEAR_ZERO);

    if (map == 1) {
      EXPECT_NEAR(pe_original,
                  system.potential(0).visit_model()->inner()->energy_map()->total_energy(),
                  10*NEAR_ZERO);
    }

    // prep to add a third particle
    Position position;
    position.set_vector({0, 1.25, 0});
    auto tsel_ghost = MakeTrialSelectParticle({{"particle_type", "0"}});
    // precompute adds argument {"ghost", "true"} to tsel
    PerturbAdd add;
    add.precompute(tsel_ghost.get(), &system);
    RandomMT19937 random;
    tsel_ghost->select(Select(), &system, &random);
    add.add(&system, tsel_ghost.get(), &random, position);
    EXPECT_EQ(system.configuration().num_particles(), 3);
    const double tri_distance = sqrt(1.25*1.25 + 1.25*1.25);
    const double peTri = 4*(pow(tri_distance, -12) - pow(tri_distance, -6));
    EXPECT_NEAR(2*pe_original + peTri, system.energy(), NEAR_ZERO);

    if (map == 1) {
      EXPECT_NEAR(2*pe_original + peTri,
                  system.potential(0).visit_model()->inner()->energy_map()->total_energy(),
                  10*NEAR_ZERO);
    }

    // revert the addition of the third particle
    add.revert(&system);
    EXPECT_EQ(system.configuration().num_particles(), 2);
    EXPECT_NEAR(pe_original, system.energy(), NEAR_ZERO);

    if (map == 1) {
      EXPECT_NEAR(pe_original,
                  system.potential(0).visit_model()->inner()->energy_map()->total_energy(),
                  10*NEAR_ZERO);
    }

    // translate a particle
    Position trajectory({0., 0., 0.});
    trajectory.set_coord(2, 1.25);
    PerturbTranslate translate;
    auto tsel = MakeTrialSelectParticle({{"particle_type", "0"}});
    tsel->select(Select(), &system, &random);
    translate.perturb(&system, tsel.get(), &random);
    EXPECT_NE(pe_original, system.energy());
    if (map == 1) {
      EXPECT_NE(pe_original,
                system.potential(0).visit_model()->inner()->energy_map()->total_energy());
    }
    translate.revert(&system);
    EXPECT_NEAR(pe_original, system.energy(), NEAR_ZERO);

    if (map == 1) {
      EXPECT_NEAR(pe_original,
                  system.potential(0).visit_model()->inner()->energy_map()->total_energy(),
                  10*NEAR_ZERO);
    }

    // remove one of the remaining two particles
    PerturbRemove remove;
    tsel->select(Select(), &system, &random);
    remove.perturb(&system, tsel.get(), &random);
    EXPECT_EQ(2, system.configuration().num_particles());

    // revert removal
    remove.revert(&system);

    EXPECT_NEAR(pe_original, system.energy(), NEAR_ZERO);
    if (map == 1) {
      EXPECT_NEAR(pe_original,
                  system.potential(0).visit_model()->inner()->energy_map()->total_energy(),
                  10*NEAR_ZERO);
    }

    // finalize removal
    remove.perturb(&system, tsel.get(), &random);
    remove.finalize(&system);
    EXPECT_EQ(1, system.configuration().num_particles());
    EXPECT_EQ(0., system.energy());

    // uncomment this to fix map reversion upon particle removal
    if (map == 1) {
      EXPECT_NEAR(0.,
                  system.potential(0).visit_model()->inner()->energy_map()->total_energy(),
                  10*NEAR_ZERO);
    }
  }
}

}  // namespace feasst

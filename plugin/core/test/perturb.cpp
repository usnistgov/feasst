#include <gtest/gtest.h>
#include <random>
#include "core/include/perturb_translate.h"
#include "core/include/perturb_transfer.h"
#include "core/include/model_lj.h"

TEST(Perturb, Revert) {
  feasst::seed_random_by_date();
  feasst::System system;
  system.default_system();
  feasst::PerturbTransfer attempt;
  feasst::Particle particle;
  particle.default_particle();
  particle.set_coordinate(1, 1.25);
  feasst::Site site = particle.site(0);
  site.set_coordinate(1, 1.25);
  particle.set_site(0, site);
  feasst::System system2(system);
  feasst::ModelLJ model;
  feasst::VisitModel visit;
  feasst::Configuration* config = system.configuration(0);
  model.compute(visit, *config, -1);
  const double peOriginal = 4*(pow(1.25, -12) - pow(1.25, -6));
  EXPECT_NEAR(peOriginal, visit.energy(), 1e-15);
  attempt.add(particle, &system);
  EXPECT_EQ(system.num_particles(), 3);
  model.compute(visit, *config, -1);
  const double triDist = sqrt(1.25*1.25 + 1.25*1.25);
  const double peTri = 4*(pow(triDist, -12) - pow(triDist, -6));
  EXPECT_NEAR(2*peOriginal + peTri, visit.energy(), 1e-15);
  attempt.revert();
  EXPECT_EQ(system.num_particles(), 2);
  model.compute(visit, *config, -1);
  EXPECT_NEAR(peOriginal, visit.energy(), 1e-15);

  model.compute(visit, *config, -1);
  EXPECT_NEAR(peOriginal, visit.energy(), 1e-15);
  feasst::Configuration* config2 = system2.configuration(0);
  //config2->select_particle(1);
  attempt.select_random_particle(0, config);
  attempt.remove_selected_particle(&system2);
  EXPECT_EQ(1, system2.num_particles());
  model.compute(visit, *config2, -1);
  EXPECT_EQ(0., visit.energy());
  attempt.revert();
  EXPECT_EQ(2, system2.num_particles());
  model.compute(visit, *config, -1);
  EXPECT_NEAR(peOriginal, visit.energy(), 1e-15);

  feasst::Position trajectory({0., 0., 0.});
  // trajectory.set_to_origin_3D();
  trajectory.set_coord(2, 1.25);
  feasst::PerturbTranslate attempt2;
  attempt2.select_random_particle(0, *config);
  attempt2.translate_selected_particle(trajectory, &system2);
  // EXPECT_EQ(2, system2.num_particles());
  // model.compute(visit, *config2, -1);
  // EXPECT_NEAR(peTri, visit.energy(), 1e-15);
}


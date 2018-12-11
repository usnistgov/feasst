#include <gtest/gtest.h>
#include "core/include/model_lj.h"
#include "core/include/model_lrc.h"
#include "core/include/visit_model.h"
#include "core/include/constants.h"
#include "core/include/file_lmp.h"
#include "core/include/file_xyz.h"
#include "core/include/physical_constants.h"

TEST(VisitModel, energy) {
  feasst::Configuration config;
  config.default_configuration();
  const double pos = 1.25;
  EXPECT_EQ(config.particle(0).position().coord(0), 0);
  EXPECT_EQ(config.particle(1).position().coord(0), pos);
  feasst::ModelLJ model;
  feasst::VisitModel visit;
  visit.loop_by_particle(config, model);
  EXPECT_NEAR(4*(pow(pos, -12) - pow(pos, -6)), visit.energy(), 1e-15);

  // check PBCs
  feasst::Position position = config.particle(1).position();
  position.set_coord(0, 3);

  feasst::Particle particle = config.particle(1);
  particle.set_position(position);
  feasst::Site site = particle.site(0);
  site.set_position(position);
  particle.set_site(0, site);
  config.select_particle(1);
  config.replace_selected_particle_position(particle);
  EXPECT_EQ(3, config.particle(0).site(0).position().size());
  EXPECT_EQ(0, config.particle(0).site(0).position().coord(0));
  EXPECT_EQ(3, config.particle(1).site(0).position().coord(0));

  EXPECT_EQ(config.particle(1).position().coord(0), 3);
  visit.loop_by_particle(config, model);
  EXPECT_NEAR(4*(pow(2, -12) - pow(2, -6)), visit.energy(), 1e-15);
  config.select_particle(0);
  //visit.loop_by_particle(config, model, 0);
  visit.energy_of_selection(config, model);
  EXPECT_NEAR(4*(pow(2, -12) - pow(2, -6)), visit.energy(), 1e-15);
}

TEST(VisitModel, reference_config) {
  feasst::Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  feasst::FileXYZ().load("../plugin/core/test/data/lj_sample_config_periodic4.xyz", &config);
  EXPECT_EQ(1, config.particle_types().num());
  EXPECT_EQ(1, config.particle_type(0).num_sites());
  feasst::ModelLJ model;
  feasst::VisitModel visit;
  visit.loop_by_particle(config, model);
  EXPECT_NEAR(-16.790321304625856, visit.energy(), 1e-15);
  feasst::ModelLRC lrc;
  visit.loop_by_particle(config, lrc);
  EXPECT_NEAR(-0.5451660014945704, visit.energy(), 1e-15);
}

TEST(VisitModel, ModelLRC) {
  feasst::Configuration config;
  config.default_configuration();
  feasst::VisitModel visit;
  feasst::ModelLRC lrc;
  visit.loop_by_particle(config, lrc);
  const double pe_lrc = (8./3.)*feasst::PI*pow(config.num_particles(), 2)/config.domain().volume()
    *((1./3.)*pow(3, -9) - pow(3, -3));
  EXPECT_NEAR(pe_lrc, visit.energy(), 1e-15);

  // test visit design pattern
  feasst::Model* model = &lrc;
  EXPECT_NEAR(pe_lrc, model->compute(visit, config, -1), 1e-15);
}

TEST(VisitModel, spce_reference_config) {
  feasst::Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  EXPECT_EQ(1, config.particle_types().num());
  EXPECT_EQ(3, config.particle_type(0).num_sites());
  for (int part = 0; part < 100; ++part) {
    config.add_particle(0);
  }
  feasst::FileXYZ().load("../plugin/core/test/data/spce_sample_config_periodic1.xyz", &config);
  EXPECT_NEAR(8000, config.domain().volume(), 1e-15);
  feasst::ModelLJ model;
  feasst::VisitModel visit;
  visit.loop_by_particle(config, model);
  EXPECT_NEAR(99538.736236886805*feasst::ideal_gas_constant/1e3, visit.energy(), 1e-15);
  feasst::ModelLRC lrc;
  visit.loop_by_particle(config, lrc);
  EXPECT_NEAR(-823.71499511652326*feasst::ideal_gas_constant/1e3, visit.energy(), 1e-13);
}


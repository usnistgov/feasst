#include <gtest/gtest.h>
#include "core/include/model_lj.h"
#include "core/include/model_lrc.h"
#include "core/include/visit_model.h"
#include "core/include/constants.h"
#include "core/include/file_lmp.h"
#include "core/include/file_xyz.h"
#include "core/include/physical_constants.h"
#include "core/include/select_list.h"

TEST(VisitModel, energy) {
  feasst::Configuration config;
  config.default_configuration();
  const double pos = 1.25;
  EXPECT_EQ(config.particle(0).position().coord(0), 0);
  EXPECT_EQ(config.particle(1).position().coord(0), pos);
  feasst::ModelLJ model;
  feasst::VisitModel visit;
  visit.compute(config, model);
  EXPECT_NEAR(4*(pow(pos, -12) - pow(pos, -6)), visit.energy(), feasst::NEAR_ZERO);

  // check PBCs
  feasst::Position position = config.particle(1).position();
  position.set_coord(0, 3);

  feasst::Particle particle = config.particle(1);
  particle.set_position(position);
  feasst::Site site = particle.site(0);
  site.set_position(position);
  particle.set_site(0, site);
  feasst::SelectList select;
  select.particle(1, config);
  config.replace_position(select, particle);
  // config.select_particle(1);
  // config.replace_selected_particle_position(particle);
  EXPECT_EQ(3, config.particle(0).site(0).position().size());
  EXPECT_EQ(0, config.particle(0).site(0).position().coord(0));
  EXPECT_EQ(3, config.particle(1).site(0).position().coord(0));

  EXPECT_EQ(config.particle(1).position().coord(0), 3);
  model.compute(visit, config, select);
  EXPECT_NEAR(4*(pow(2, -12) - pow(2, -6)), visit.energy(), feasst::NEAR_ZERO);
  config.select_particle(0);
  select.particle(0, config);
  model.compute(visit, config, select);
  EXPECT_NEAR(4*(pow(2, -12) - pow(2, -6)), visit.energy(), feasst::NEAR_ZERO);
}

TEST(VisitModel, reference_config) {
  feasst::Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  feasst::FileXYZ().load("../plugin/core/test/data/lj_sample_config_periodic4.xyz", &config);
  EXPECT_EQ(1, config.particle_types().num());
  EXPECT_EQ(1, config.particle_type(0).num_sites());
  feasst::ModelLJ model;
  feasst::VisitModel visit;
  visit.compute(config, model);
  EXPECT_NEAR(-16.790321304625856, visit.energy(), feasst::NEAR_ZERO);
  feasst::ModelLRC lrc;
  visit.compute(config, lrc);
  EXPECT_NEAR(-0.5451660014945704, visit.energy(), feasst::NEAR_ZERO);
  visit.check_energy(config, model);
}

TEST(VisitModel, ModelLRC) {
  feasst::Configuration config;
  config.default_configuration();
  feasst::VisitModel visit;
  feasst::ModelLRC lrc;
  visit.compute(config, lrc);
  const double pe_lrc = (8./3.)*feasst::PI*pow(config.num_particles(), 2)/config.domain().volume()
    *((1./3.)*pow(3, -9) - pow(3, -3));
  EXPECT_NEAR(pe_lrc, visit.energy(), feasst::NEAR_ZERO);

  // test visit design pattern
  feasst::Model* model = &lrc;
  EXPECT_NEAR(pe_lrc, model->compute(visit, config), feasst::NEAR_ZERO);
}

TEST(VisitModel, spce_reference_config) {
  feasst::seed_random_by_date();
  feasst::Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  EXPECT_EQ(1, config.particle_types().num());
  EXPECT_EQ(3, config.particle_type(0).num_sites());
  for (int part = 0; part < 100; ++part) {
    config.add_particle(0);
  }
  feasst::FileXYZ().load("../plugin/core/test/data/spce_sample_config_periodic1.xyz", &config);
  EXPECT_NEAR(8000, config.domain().volume(), feasst::NEAR_ZERO);
  feasst::ModelLJ model;
  feasst::VisitModel visit;
  visit.compute(config, model);
  const double pe_lj = 99538.736236886805;
  EXPECT_NEAR(pe_lj*feasst::ideal_gas_constant/1e3, visit.energy(), feasst::NEAR_ZERO);
  feasst::ModelLRC lrc;
  visit.compute(config, lrc);
  const double pe_lrc = -823.71499511652326;
  EXPECT_NEAR(pe_lrc*feasst::ideal_gas_constant/1e3, visit.energy(), 1e-13);

  // test adding/deleting particles, resulting in a ghost
  feasst::SelectList select;
  select.random_particle(config);
  visit.compute(config, model, select);
  const double pe_previous = visit.energy();
  const double x1_previous = select.particle(config).site(0).position().coord(1);
  config.add_particle(0);
  feasst::SelectList new_part;
  new_part.last_particle_added(&config);
  config.replace_position(new_part, select.particle(config));
  config.remove_particle(select);
  visit.compute(config, model);
  EXPECT_NEAR(pe_lj*feasst::ideal_gas_constant/1e3, visit.energy(), 1e-12);
  visit.compute(config, lrc);
  EXPECT_NEAR(pe_lrc*feasst::ideal_gas_constant/1e3, visit.energy(), 1e-13);
  EXPECT_EQ(101, config.particles().num()); // includes one ghost particle
  EXPECT_EQ(100, config.selection_of_all().num_particles());
  config.check_size();
  visit.compute(config, model, new_part);
  EXPECT_NEAR(pe_previous, visit.energy(), feasst::NEAR_ZERO);
  const double x2_previous = new_part.particle(config).site(0).position().coord(1);
  EXPECT_EQ(x1_previous, x2_previous);
  visit.check_energy(config, model);
}

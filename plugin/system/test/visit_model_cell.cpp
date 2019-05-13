#include "utils/test/utils.h"
#include "configuration/test/configuration_test.h"
#include "system/include/visit_model_cell.h"
#include "configuration/include/file_xyz.h"
#include "system/include/model_lj.h"
#include "math/include/constants.h"
#include "system/include/physical_constants.h"
#include "system/include/select_list.h"
#include "math/include/utils_math.h"
#include "monte_carlo/include/perturb_configs.h"

namespace feasst {

/// Use a 5 angstrom cut off
TEST(VisitModelCell, simple_lj) {
  seed_random_by_date();
  Configuration config;
  config.set_domain(Domain().set_cubic(15));
  config.add_particle_type("../forcefield/data.lj");
  config.add_particle_of_type(0);
  SelectList select;
  config.add_particle_of_type(0);
  select.last_particle_added(&config);
  auto select2 = select;
  config.add_particle_of_type(0);
  select2.last_particle_added(&config);
  config.remove_particle(select2);
  feasst::Position pos;
  pos.set_vector({2, 0, 0});
  config.displace_particle(select, pos);
  EXPECT_EQ(0, config.particle(0).site(0).position().coord(0));
  EXPECT_EQ(2, config.particle(1).site(0).position().coord(0));
  config.init_cells(3);
  const Cells& cells = config.domain().cells(0);
  EXPECT_EQ(5*5*5, cells.num_total());
  const int center = round(5.*5.*5./2. - 0.5);
  EXPECT_EQ(config.particle(0).site(0).property("cell0"), center);
  EXPECT_EQ(config.particle(1).site(0).property("cell0"), center + 1);
  EXPECT_EQ(cells.particles()[center].num_sites(), 1);
  EXPECT_EQ(cells.particles()[center].particle_index(0), 0);
  EXPECT_EQ(cells.particles()[center].site_index(0, 0), 0);
  EXPECT_EQ(cells.particles()[center + 1].num_sites(), 1);
  EXPECT_EQ(cells.particles()[center + 1].particle_index(0), 1);
  EXPECT_EQ(cells.particles()[center + 1].site_index(0, 0), 0);
  config.check();
  ModelLJ model;
  VisitModelCell cell_visit;
  VisitModel visit;
  model.compute(&config, &visit);
  model.compute(&config, &cell_visit);
  double r2 = 4;
  double pe_lj = 4*(pow(r2, -6) - pow(r2, -3));
  EXPECT_NEAR(visit.energy(), pe_lj, NEAR_ZERO);
  EXPECT_NEAR(visit.energy(), cell_visit.energy(), NEAR_ZERO);

  visit.check_energy(model, &config);
  cell_visit.check_energy(model, &config);

  DEBUG("Move to the same cell");
  pos.set_vector({-1, 1, 1});
  config.displace_particle(select, pos);
  EXPECT_EQ(config.particle(0).site(0).property("cell0"), center);
  EXPECT_EQ(config.particle(1).site(0).property("cell0"), center);
  EXPECT_EQ(cells.particles()[center].num_sites(), 2);
  EXPECT_EQ(cells.particles()[center].num_particles(), 2);
  EXPECT_EQ(cells.particles()[center].particle_index(0), 0);
  EXPECT_EQ(cells.particles()[center].particle_index(1), 1);
  EXPECT_EQ(cells.particles()[center].site_index(0, 0), 0);
  EXPECT_EQ(cells.particles()[center].site_index(1, 0), 0);
  model.compute(&config, &visit);
  model.compute(&config, &cell_visit);
  r2 = 3;
  pe_lj = 4*(pow(r2, -6) - pow(r2, -3));
  EXPECT_NEAR(visit.energy(), pe_lj, NEAR_ZERO);
  EXPECT_NEAR(visit.energy(), cell_visit.energy(), NEAR_ZERO);

  /// test energy of selection
  model.compute(select, &config, &visit);
  model.compute(select, &config, &cell_visit);
  EXPECT_NEAR(visit.energy(), pe_lj, NEAR_ZERO);
  EXPECT_NEAR(visit.energy(), cell_visit.energy(), 5e-15);

  visit.check_energy(model, &config);
  cell_visit.check_energy(model, &config);
}

/// Use a 5 angstrom cut off
TEST(VisitModelCell, lj_reference_config) {
  seed_random_by_date();
  Configuration config = lj_sample();
  const int rcut = 2;
  for (int site_index = 0; site_index < config.num_site_types(); ++site_index) {
    config.set_model_param("cutoff", site_index, rcut);
  }
  config.add_particle_of_type(0);
  SelectList select;
  select.last_particle_added(&config);
  config.remove_particle(select);
  config.init_cells(rcut);
  config.check();
  ModelLJ model;
  VisitModelCell cell_visit;
  VisitModel visit;
  model.compute(&config, &visit);
  model.compute(&config, &cell_visit);
  EXPECT_NEAR(visit.energy(), cell_visit.energy(), 5e-12);
  EXPECT_NEAR(-15.076312312129405, visit.energy(), NEAR_ZERO);

  /// test energy of selection
  select.random_particle_of_type(0, &config);
  model.compute(select, &config, &visit);
  model.compute(select, &config, &cell_visit);
  EXPECT_NEAR(visit.energy(), cell_visit.energy(), 5e-15);

  auto cell_visit2 = test_serialize<VisitModelCell, VisitModel>(cell_visit);
  model.compute(select, &config, cell_visit2.get());
  EXPECT_NEAR(visit.energy(), cell_visit2->energy(), 5e-15);

  visit.check_energy(model, &config);
  cell_visit.check_energy(model, &config);
}

/// Use a 5 angstrom cut off
TEST(VisitModelCell, spce_reference_config) {
  seed_random_by_date();
  // seed_random();
  Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  const int rcut = 5;
  for (int site_index = 0; site_index < config.num_site_types(); ++site_index) {
    config.set_model_param("cutoff", site_index, rcut);
  }
  try {
    auto config2 = config;
    config2.init_cells(rcut);
    CATCH_PHRASE("cannot define cells before domain side");
  }
  for (int part = 0; part < 100; ++part) {
    config.add_particle_of_type(0);
  }
  FileXYZ().load("../plugin/system/test/data/spce_sample_config_periodic1.xyz", &config);
  config.add_particle_of_type(0);
  SelectList select;
  select.last_particle_added(&config);
  config.remove_particle(select);
  config.init_cells(rcut);
  config.check();
  ModelLJ model;
  VisitModelCell cell_visit;
  VisitModel visit;
  model.compute(&config, &visit);
  model.compute(&config, &cell_visit);
  EXPECT_NEAR(visit.energy(), cell_visit.energy(), 5e-12);
  EXPECT_NEAR(896.85497602741475, visit.energy(), NEAR_ZERO);

  /// test energy of selection
  select.random_particle_of_type(0, &config);
  model.compute(select, &config, &visit);
  model.compute(select, &config, &cell_visit);
  EXPECT_NEAR(visit.energy(), cell_visit.energy(), 5e-14);

  visit.check_energy(model, &config);
  cell_visit.check_energy(model, &config);
}

// add individual particles until reaching the point where the visitors are
// inconsistent
TEST(VisitModelCell, spce_reference_config_buildup) {
  seed_random_by_date();
  // seed_random();
  Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  const int rcut = 5;
  for (int site_index = 0; site_index < config.num_site_types(); ++site_index) {
    config.set_model_param("cutoff", site_index, rcut);
  }
  for (int part = 0; part < 100; ++part) {
    config.add_particle_of_type(0);
  }
  FileXYZ().load("../plugin/system/test/data/spce_sample_config_periodic1.xyz", &config);
  config.init_cells(rcut);
  ModelLJ model;
  VisitModelCell cell_visit;
  VisitModel visit;

  System sys;
  sys.add(config);
  sys.add(config);
  Configuration * config1 = sys.get_configuration(0);
  Configuration * config2 = sys.get_configuration(1);

  // remove particles in config2
  EXPECT_EQ(100, config1->num_particles());
  EXPECT_EQ(100, config2->num_particles());
  config2->remove_particles(config2->group_select(0));
  EXPECT_EQ(0, config2->num_particles());

  PerturbConfigs perturb;
  int transfers = 0;
  while (config1->num_particles() > 90) {
    perturb.transfer_particle(0, &sys, 0, 1);
    ++transfers;
    DEBUG("transfers " << transfers);
    DEBUG("cell list " << config2->domain().cells(0).str());
    EXPECT_EQ(100 - transfers, config1->num_particles());
    EXPECT_EQ(transfers, config2->num_particles());
    cell_visit.check_energy(model, config2);
  }
}

}  // namespace feasst


#include <gtest/gtest.h>
#include "core/include/visit_model_cell.h"
#include "core/include/file_xyz.h"
#include "core/include/model_lj.h"
#include "core/include/constants.h"
#include "core/include/physical_constants.h"
#include "core/include/select_list.h"
#include "core/include/utils_math.h"

/// Use a 5 angstrom cut off
TEST(VisitModelCell, simple_lj) {
  feasst::seed_random_by_date();
  feasst::Configuration config;
  config.set_domain(feasst::Domain().set_cubic(15));
  config.add_particle_type("../forcefield/data.lj");
  config.add_particle(0);
  feasst::SelectList select;
  config.add_particle(0);
  select.select_sites(config, 1, {0});
  feasst::Position pos;
  pos.set_vector({2, 0, 0});
  config.displace_particle(select, pos);
  const feasst::Site& site0 = config.particle(0).site(0);
  EXPECT_EQ(0, site0.position().coord(0));
  const feasst::Site& site1 = config.particle(1).site(0);
  EXPECT_EQ(2, site1.position().coord(0));
  config.init_cells(3);
  const feasst::Cells& cells = config.domain().cells(0);
  EXPECT_EQ(5*5*5, cells.num_total());
  const int center = feasst::round(5.*5.*5./2. - 0.5);
  EXPECT_EQ(site0.property("cell0"), center);
  EXPECT_EQ(site1.property("cell0"), center + 1);
  EXPECT_EQ(cells.particles()[center].num_sites(), 1);
  EXPECT_EQ(cells.particles()[center].particle_index(0), 0);
  EXPECT_EQ(cells.particles()[center].site_index(0, 0), 0);
  EXPECT_EQ(cells.particles()[center + 1].num_sites(), 1);
  EXPECT_EQ(cells.particles()[center + 1].particle_index(0), 1);
  EXPECT_EQ(cells.particles()[center + 1].site_index(0, 0), 0);
  config.check_size();
  feasst::ModelLJ model;
  feasst::VisitModelCell cell_visit;
  feasst::VisitModel visit;
  model.compute(visit, config);
  model.compute(cell_visit, config);
  double r2 = 4;
  double pe_lj = 4*(pow(r2, -6) - pow(r2, -3));
  EXPECT_NEAR(visit.energy(), pe_lj, feasst::NEAR_ZERO);
  EXPECT_NEAR(visit.energy(), cell_visit.energy(), feasst::NEAR_ZERO);

  DEBUG("Move to the same cell");
  pos.set_vector({-1, 1, 1});
  config.displace_particle(select, pos);
  EXPECT_EQ(site0.property("cell0"), center);
  EXPECT_EQ(site1.property("cell0"), center);
  EXPECT_EQ(cells.particles()[center].num_sites(), 2);
  EXPECT_EQ(cells.particles()[center].num_particles(), 2);
  EXPECT_EQ(cells.particles()[center].particle_index(0), 0);
  EXPECT_EQ(cells.particles()[center].particle_index(1), 1);
  EXPECT_EQ(cells.particles()[center].site_index(0, 0), 0);
  EXPECT_EQ(cells.particles()[center].site_index(1, 0), 0);
  model.compute(visit, config);
  model.compute(cell_visit, config);
  r2 = 3;
  pe_lj = 4*(pow(r2, -6) - pow(r2, -3));
  EXPECT_NEAR(visit.energy(), pe_lj, feasst::NEAR_ZERO);
  EXPECT_NEAR(visit.energy(), cell_visit.energy(), feasst::NEAR_ZERO);
}

/// Use a 5 angstrom cut off
TEST(VisitModelCell, spce_reference_config) {
  feasst::seed_random_by_date();
  feasst::Configuration config;
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
    config.add_particle(0);
  }
  feasst::FileXYZ().load("../plugin/core/test/data/spce_sample_config_periodic1.xyz", &config);
  config.init_cells(rcut);
  config.check_size();
  feasst::ModelLJ model;
  feasst::VisitModelCell cell_visit;
  feasst::VisitModel visit;
  model.compute(visit, config);
  model.compute(cell_visit, config);
  EXPECT_NEAR(visit.energy(), cell_visit.energy(), 5e-12);
  //const double pe_lj = 99538.736236886805;
  //EXPECT_NEAR(pe_lj*feasst::ideal_gas_constant/1e3, visit.energy(), feasst::NEAR_ZERO);
}

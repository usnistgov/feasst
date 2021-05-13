#include "utils/test/utils.h"
#include "math/include/constants.h"
#include "math/include/random_mt19937.h"
#include "math/include/utils_math.h"
#include "system/include/visit_model_cell.h"
#include "system/include/lennard_jones.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/domain.h"
#include "configuration/include/select.h"

namespace feasst {

TEST(VisitModelCell, cells) {
  Configuration config(MakeDomain({{"cubic_box_length", "7"}}),
    {{"particle_type0", "../forcefield/data.spce"}});
  config.add(MakeGroup({{"site_type", "0"}}));
  config.add_particle_of_type(0);

  auto visit = MakeVisitModelCell({{"min_length", "1.4"}});
  visit->precompute(&config);

  EXPECT_EQ(visit->cells().num_total(), 5*5*5);
//  int cell0 = round(7*7*7/2. - 0.5);
  int cell1 = round(5*5*5/2. - 0.5);
  Site site = config.particle(0).site(0);
//  EXPECT_EQ(cell0, site.cell(0));
  EXPECT_EQ(cell1, site.cell(0));
  std::vector<int> indices = {0};
//  EXPECT_EQ(config.domain().cells(0).particles()[0].num_particles(), 0);
  EXPECT_EQ(visit->cells().particles()[0].num_particles(), 0);
//  EXPECT_EQ(config.domain().cells(0).particles()[cell0].num_particles(), 1);
  EXPECT_EQ(visit->cells().particles()[cell1].num_particles(), 1);
//  EXPECT_EQ(config.particle(0).site(0).num_cells(), 2);
  EXPECT_EQ(config.particle(0).site(1).num_cells(), 1);
  Position trajectory({-3.49, -3.49, -3.49});

  DEBUG("displacing particles");
  Select select(0, config.select_particle(0));
  //select.particle(0, config);
  config.displace_particles(select, trajectory);
  visit->finalize(select);
  EXPECT_EQ(0, config.particle(0).site(0).cell(0));
  site = config.particle(0).site(0);
  EXPECT_EQ(0, site.cell(0));
//  EXPECT_EQ(0, site.cell(1));
  DEBUG("pos " << config.particle(0).site(2).position().coord(0));
  //EXPECT_EQ(1, config.particle(0).site(1).cell(0));
  //EXPECT_EQ(6, config.particle(0).site(2).cell(0));
  EXPECT_EQ(config.particle(0).site(0).num_cells(), 1);
  EXPECT_EQ(config.particle(0).site(1).num_cells(), 1);
  EXPECT_EQ(config.particle(0).site(2).num_cells(), 1);
  EXPECT_EQ(visit->cells().particles()[cell1].num_particles(), 0);
//  EXPECT_EQ(config.domain().cells(1).particles()[cell1].num_particles(), 0);
  EXPECT_NE(cell1, site.cell(0));

  // serialize
//  Configuration config2 = test_serialize(config);
//  EXPECT_EQ(6, config2.particle(0).site(2).cell(0));
//  EXPECT_EQ(config2.domain().cells(1).particles()[cell1].num_particles(), 0);
//
//  config.remove_particle(select);
  visit->check();
}

/// Use a 5 angstrom cut off
TEST(VisitModelCell, simple_lj) {
  Configuration config(MakeDomain({{"cubic_box_length", "15"}}),
    {{"particle_type", "../forcefield/data.lj"}});
  config.add_particle_of_type(0);
  config.add_particle_of_type(0);
  Select select(config.newest_particle_index(), config.newest_particle());
  config.add_particle_of_type(0);
  Select select2(config.newest_particle_index(), config.newest_particle());
  config.remove_particle(select2);
  feasst::Position pos;
  pos.set_vector({2, 0, 0});
  config.displace_particle(select, pos);
  EXPECT_EQ(0, config.particle(0).site(0).position().coord(0));
  EXPECT_EQ(2, config.particle(1).site(0).position().coord(0));
  LennardJones model;
  auto cell_visit = MakeVisitModelCell({{"min_length", "3"}});
  cell_visit->precompute(&config);
//  config.init_cells(3);
  const Cells& cells = cell_visit->cells();
  EXPECT_EQ(5*5*5, cells.num_total());
  const int center = round(5.*5.*5./2. - 0.5);
  EXPECT_EQ(config.particle(0).site(0).cell(0), center);
  EXPECT_EQ(config.particle(1).site(0).cell(0), center + 1);
  EXPECT_EQ(cells.particles()[center].num_sites(), 1);
  EXPECT_EQ(cells.particles()[center].particle_index(0), 0);
  EXPECT_EQ(cells.particles()[center].site_index(0, 0), 0);
  EXPECT_EQ(cells.particles()[center + 1].num_sites(), 1);
  EXPECT_EQ(cells.particles()[center + 1].particle_index(0), 1);
  EXPECT_EQ(cells.particles()[center + 1].site_index(0, 0), 0);
  config.check();
  cell_visit->check();
  VisitModel visit;
  visit.precompute(&config);
  model.compute(&config, &visit);
  model.compute(&config, cell_visit.get());
  double r2 = 4;
  double pe_lj = 4*(pow(r2, -6) - pow(r2, -3));
  EXPECT_NEAR(visit.energy(), pe_lj, NEAR_ZERO);
  EXPECT_NEAR(visit.energy(), cell_visit->energy(), NEAR_ZERO);

  visit.check_energy(&model, &config);
  cell_visit->check_energy(&model, &config);

  DEBUG("Move to the same cell");
  pos.set_vector({-1, 1, 1});
  config.displace_particle(select, pos);
  cell_visit->finalize(select);
  EXPECT_EQ(config.particle(0).site(0).cell(0), center);
  EXPECT_EQ(config.particle(1).site(0).cell(0), center);
  EXPECT_EQ(cells.particles()[center].num_sites(), 2);
  EXPECT_EQ(cells.particles()[center].num_particles(), 2);
  EXPECT_EQ(cells.particles()[center].particle_index(0), 0);
  EXPECT_EQ(cells.particles()[center].particle_index(1), 1);
  EXPECT_EQ(cells.particles()[center].site_index(0, 0), 0);
  EXPECT_EQ(cells.particles()[center].site_index(1, 0), 0);
  model.compute(&config, &visit);
  model.compute(&config, cell_visit.get());
  r2 = 3;
  pe_lj = 4*(pow(r2, -6) - pow(r2, -3));
  EXPECT_NEAR(visit.energy(), pe_lj, NEAR_ZERO);
  EXPECT_NEAR(visit.energy(), cell_visit->energy(), NEAR_ZERO);

  /// test energy of selection
  model.compute(select, &config, &visit);
  model.compute(select, &config, cell_visit.get());
  EXPECT_NEAR(visit.energy(), pe_lj, NEAR_ZERO);
  EXPECT_NEAR(visit.energy(), cell_visit->energy(), 5e-15);

  visit.check_energy(&model, &config);
  cell_visit->check_energy(&model, &config);
}

}  // namespace feasst

#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "system/include/visit_model_cell.h"
#include "configuration/include/file_xyz.h"
#include "system/include/lennard_jones.h"
#include "math/include/constants.h"
#include "configuration/include/select.h"
#include "math/include/utils_math.h"

namespace feasst {

/// Use a 5 angstrom cut off
TEST(VisitModelCell, simple_lj) {
  Configuration config(MakeDomain({{"cubic_box_length", "15"}, {"init_cells", "3"}}),
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
//  config.init_cells(3);
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
  LennardJones model;
  VisitModelCell cell_visit;
  cell_visit.precompute(&config);
  VisitModel visit;
  visit.precompute(&config);
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

}  // namespace feasst

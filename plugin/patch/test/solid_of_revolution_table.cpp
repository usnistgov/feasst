#include "utils/test/utils.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/domain.h"
#include "system/include/system.h"
#include "system/include/potential.h"
#include "system/include/visit_model.h"
#include "models/include/square_well.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "patch/include/solid_of_revolution_table.h"

namespace feasst {

TEST(SolidOfRevolutionTable, serialize) {
  SolidOfRevolutionTable patch;
  SolidOfRevolutionTable patch2 = test_serialize(patch);
}

TEST(SolidOfRevolutionTable, test_two) {
  System system;
  system.add(MakeConfiguration({{"cubic_side_length", "10"},
    {"particle_type", "../plugin/patch/particle/one_patch.fstprt"},
    {"add_particles_of_type0", "2"},
    {"group0", "centers"}, {"centers_site_type0", "0"}}));

  system.add(MakePotential(MakeSquareWell(),
                       MakeVisitModel(std::make_shared<SolidOfRevolutionTable>(argtype({{"table_file", "../plugin/patch/test/data/tablek5l2.0d1.txt"}}))),
                       //MakeVisitModel(MakeSolidOfRevolutionTable({{"table_file", "../plugin/patch/test/data/tablek5l2.0d1.txt"}})),
                       {{"group", "centers"}}));  // optimization: loop centers
  system.precompute();

  PerturbAnywhere anywhere;
  TrialSelectParticle tsel;
  tsel.select_particle(1, system.configuration());
  auto traj = Position().set_vector({1.001, 0., 0.});
  auto file = MakeFileXYZ({{"append", "true"}});
  file->write_for_vmd("tmp/sor.xyz", system.configuration());

  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sor.xyz", system.configuration());
  EXPECT_NEAR(0, system.energy(), NEAR_ZERO);

  traj = Position().set_vector({0.999, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sor.xyz", system.configuration());
  EXPECT_NEAR(NEAR_INFINITY, system.energy(), 1e100);

  // rotate 45 degrees
  PerturbRotate rotate;
  auto axis = Position().set_vector({0., 1., 0.});
  rotate.move(
    traj,
    RotationMatrix().axis_angle(axis, 45.),
    &system,
    &tsel
  );
//  file->write("tmp/sor.xyz", system.configuration());
//  EXPECT_NEAR(-1, system.energy(), NEAR_ZERO);

  // orientation: |/ with 45 degree angle contacts at L/2+D/sqrt(2)
  double sep = 1+1./std::sqrt(2);
  traj = Position().set_vector({sep + 1e-4, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sor.xyz", system.configuration());
  EXPECT_NEAR(0, system.energy(), NEAR_ZERO);

  traj = Position().set_vector({sep - 1e-4, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sor.xyz", system.configuration());
  EXPECT_NEAR(NEAR_INFINITY, system.energy(), 1e100);

  // rotate another 45 degrees
  rotate.move(
    traj,
    RotationMatrix().axis_angle(axis, 45.),
    &system,
    &tsel
  );
  file->write("tmp/sor.xyz", system.configuration());

  traj = Position().set_vector({2.001, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sor.xyz", system.configuration());
  EXPECT_NEAR(0, system.energy(), NEAR_ZERO);

  traj = Position().set_vector({1.999, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sor.xyz", system.configuration());
  EXPECT_NEAR(NEAR_INFINITY, system.energy(), 1e100);
}

} // namespace feasst


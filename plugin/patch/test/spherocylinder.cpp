#include "utils/test/utils.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/domain.h"
#include "configuration/include/group.h"
#include "system/include/visit_model.h"
#include "system/include/system.h"
#include "system/include/potential.h"
#include "models/include/square_well.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "patch/include/spherocylinder.h"

namespace feasst {

TEST(Spherocylinder, serialize) {
  Spherocylinder patch;
  Spherocylinder patch2 = test_serialize(patch);
}

TEST(Spherocylinder, patch_one) {
  Configuration config;
  config.add_particle_type("../plugin/patch/particle/spherocylinder.fstprt");
  config.set_model_param("cutoff", 0, 3.);
  config.set_model_param("cutoff", 1, 3.);
  FileXYZ().load("../plugin/patch/test/data/spherocylinder.xyz", &config);
  SquareWell model;
  model.precompute(config.model_params());
  VisitModel visit;
  auto patch = MakeSpherocylinder();
  visit.set_inner(patch);
  visit.precompute(&config);
  //patch->set_patch_angle(1, 90.);
  config.add(MakeGroup({{"site_type", "0"}}));
  visit.compute(&model, &config, 1);
  EXPECT_NEAR(-1., visit.energy(), NEAR_ZERO);

  auto visit2 = test_serialize<VisitModel, VisitModel>(visit);
  visit2->compute(&model, &config, 1);
  EXPECT_NEAR(-1., visit2->energy(), NEAR_ZERO);
}

TEST(Spherocylinder, patch_one_2body) {
  System system;
  system.add(MakeConfiguration({{"cubic_side_length", "10"},
    {"particle_type", "../plugin/patch/particle/spherocylinder.fstprt"},
    {"add_particles_of_type0", "2"},
    {"group0", "centers"}, {"centers_site_type0", "0"}}));

  system.add(MakePotential(MakeSquareWell(),
                       MakeVisitModel(MakeSpherocylinder()),
                       {{"group", "centers"}}));  // optimization: loop centers
  system.precompute();

  PerturbAnywhere anywhere;
  TrialSelectParticle tsel;
  tsel.select_particle(1, system.configuration());
  auto traj = Position().set_vector({3.501, 0., 0.});
  auto file = MakeFileXYZ({{"append", "true"}});
  file->write_for_vmd("tmp/sphc.xyz", system.configuration());

  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc.xyz", system.configuration());
  EXPECT_NEAR(0, system.energy(), NEAR_ZERO);

  traj = Position().set_vector({3.499, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc.xyz", system.configuration());
  EXPECT_NEAR(-1, system.energy(), NEAR_ZERO);

  traj = Position().set_vector({3.001, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc.xyz", system.configuration());
  EXPECT_NEAR(-1, system.energy(), NEAR_ZERO);

  traj = Position().set_vector({2.999, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc.xyz", system.configuration());
  EXPECT_NEAR(NEAR_INFINITY, system.energy(), 1e100);

  // rotate 45 degrees
  PerturbRotate rotate;
  auto axis = Position().set_vector({0., 0., 1.});
  rotate.move(
    traj,
    RotationMatrix().axis_angle(axis, 45.),
    &system,
    &tsel
  );
//  file->write("tmp/sphc.xyz", system.configuration());
//  EXPECT_NEAR(-1, system.energy(), NEAR_ZERO);

  // orientation: \| with 45 degree angle contacts at L/2+D/sqrt(2)
  double sep = 1+1./std::sqrt(2);
  traj = Position().set_vector({0., sep + 1e-4, 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc2.xyz", system.configuration());
  EXPECT_NEAR(-1, system.energy(), NEAR_ZERO);

  traj = Position().set_vector({0., sep - 1e-4, 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc2.xyz", system.configuration());
  EXPECT_NEAR(NEAR_INFINITY, system.energy(), 1e100);

  // orientation: -/ with 45 degree angle contacts at L/2+D*sqrt(2)
  sep = 1+std::sqrt(2);
  traj = Position().set_vector({sep + 1e-4, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc2.xyz", system.configuration());
  EXPECT_NEAR(-1, system.energy(), NEAR_ZERO);

  traj = Position().set_vector({sep - 1e-4, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc2.xyz", system.configuration());
  EXPECT_NEAR(NEAR_INFINITY, system.energy(), 1e100);

  traj = Position().set_vector({2.999, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc.xyz", system.configuration());
  EXPECT_NEAR(-1, system.energy(), NEAR_ZERO);

  // rotate another 45 degrees
  rotate.move(
    traj,
    RotationMatrix().axis_angle(axis, 45.),
    &system,
    &tsel
  );
  file->write("tmp/sphc.xyz", system.configuration());
  EXPECT_NEAR(0, system.energy(), NEAR_ZERO);

  traj = Position().set_vector({2.501, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc.xyz", system.configuration());
  EXPECT_NEAR(0, system.energy(), NEAR_ZERO);

  traj = Position().set_vector({2.499, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc.xyz", system.configuration());
  EXPECT_NEAR(-1, system.energy(), NEAR_ZERO);

  traj = Position().set_vector({2.001, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc.xyz", system.configuration());
  EXPECT_NEAR(-1, system.energy(), NEAR_ZERO);

  traj = Position().set_vector({1.999, 0., 0.});
  anywhere.set_position(traj, &system, &tsel);
  file->write("tmp/sphc.xyz", system.configuration());
  EXPECT_NEAR(NEAR_INFINITY, system.energy(), 1e100);
}

} // namespace feasst


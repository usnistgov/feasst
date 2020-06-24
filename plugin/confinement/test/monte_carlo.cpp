#include "utils/test/utils.h"
#include "utils/include/checkpoint.h"
#include "utils/include/progress_report.h"
#include "math/include/table.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "configuration/include/file_xyz.h"
#include "system/include/lennard_jones.h"
#include "shape/include/sphere.h"
#include "shape/include/slab.h"
#include "shape/include/slab_sine.h"
#include "shape/include/shape_union.h"
#include "shape/include/shape_intersect.h"
#include "shape/include/cylinder.h"
#include "shape/include/half_space_sine.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/check_energy_and_tune.h"
#include "confinement/include/model_hard_shape.h"
#include "confinement/include/model_table_cartesian.h"

namespace feasst {

TEST(MonteCarlo, ShapeUnion) {
  MonteCarlo mc;
  mc.add(Configuration(MakeDomain({{"cubic_box_length", "8"}}),
    {{"particle_type", "../forcefield/data.lj"}}));
  mc.add(Potential(MakeLennardJones()));
  mc.add(Potential(MakeModelHardShape(MakeShapeUnion(
    MakeSphere(
      {{"radius", "2"}},
      Position({{"x", "0"}, {"y", "0"}, {"z", "0"}})),
    MakeSlab({{"dimension", "2"}, {"bound0", "-1"}, {"bound1", "1"}})))));
  mc.add(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "2."}}));
  SeekNumParticles(10).with_trial_add().run(&mc);
  const int steps_per = 1e0;
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)},
                          {"file_name", "tmp/confine"}}));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(1e3);
}

std::shared_ptr<Shape> porous_network() {
  auto network = MakeShapeUnion(
    MakeSphere(
      {{"radius", "5"}},
      Position({{"x", "0"}, {"y", "0"}, {"z", "0"}})),
    MakeCylinder({{"radius", "2"}},
      Position({{"x", "0"}, {"y", "0"}, {"z", "0"}}),
      Position({{"x", "0"}, {"y", "0"}, {"z", "1"}})));
  network = MakeShapeUnion(network,
    MakeCylinder({{"radius", "2"}},
      Position({{"x", "0"}, {"y", "0"}, {"z", "0"}}),
      Position({{"x", "0"}, {"y", "1"}, {"z", "0"}})));
  network = MakeShapeUnion(network,
    MakeCylinder({{"radius", "2"}},
      Position({{"x", "0"}, {"y", "0"}, {"z", "0"}}),
      Position({{"x", "1"}, {"y", "0"}, {"z", "0"}})));
  return network;
}

TEST(MonteCarlo, ShapeUnion_LONG) {
  MonteCarlo mc;
  mc.add(Configuration(MakeDomain({{"cubic_box_length", "20"}}),
    {{"particle_type", "../forcefield/data.lj"}}));
  mc.add(Potential(MakeLennardJones()));
  mc.add(Potential(MakeModelHardShape(porous_network())));
  mc.add(MakeMetropolis({{"beta", "1.5"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "2."}}));
  SeekNumParticles(500).with_trial_add().run(&mc);
  const int steps_per = 1e5;
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)},
                          {"file_name", "tmp/confine"}}));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(1e6);
}

TEST(MonteCarlo, ShapeTable_LONG) {
  MonteCarlo mc;
  auto domain = MakeDomain({{"cubic_box_length", "20"}});
  mc.add(Configuration(domain, {{"particle_type", "../forcefield/data.lj"}}));
  mc.add(Potential(MakeLennardJones()));
  auto pore = porous_network();
  mc.add(Potential(MakeModelHardShape(pore)));
  const bool read_table = false;
  //const bool read_table = true;
  std::shared_ptr<ModelTableCart3FoldSym> hamaker;
  if (read_table) {
    auto table = MakeTable3D();
    MakeCheckpoint({{"file_name", "tmp/table"}})->read(table.get());
    hamaker = MakeModelTableCart3FoldSym(table);
  } else {
    hamaker = MakeModelTableCart3FoldSym(MakeTable3D({
      {"num0", "101"},
      {"num1", "101"},
      {"num2", "101"},
      {"default_value", "0."}}));
    #ifdef _OPENMP
    hamaker->compute_table_omp(
    #elif // _OPENMP
    hamaker->compute_table(
    #endif // _OPENMP
      pore.get(), domain.get(), mc.get_random(), {
      {"alpha", "6"},
      {"epsilon", "-1"},
      {"max_radius", "10"},
      {"num_shells", "10"},
      {"points_per_shell", "1"}});
    MakeCheckpoint({{"file_name", "tmp/table"}})->write(hamaker->table());
  }
  mc.add(Potential(hamaker));

//  EXPECT_NEAR(table2->linear_interpolation(0, 0, 0), 0., NEAR_ZERO);
//  mc.add(Potential(MakeModelTableCart3FoldSym(table2)));
//  //mc.add(Potential(MakeModelTableCart3FoldSym(table)));
//  mc.add(MakeMetropolis({{"beta", "1.5"}, {"chemical_potential", "1."}}));
//  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.1"}}));
//  SeekNumParticles(500).with_trial_add().run(&mc);
//  const int steps_per = 1e4;
//  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)},
//                          {"file_name", "tmp/confine"}}));
//  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)},
//                          {"tolerance", str(1e-4)}}));
//  FileXYZ().write("hi.xyz", mc.configuration());
//  MakeCheckpoint({{"file_name", "tmp/mc_table"}})->write(mc);
//  MonteCarlo mc2 = test_serialize(mc);
//  mc2.attempt(1e5);
}

TEST(MonteCarlo, SineSlab) {
  MonteCarlo mc;
  // mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.add(Configuration(MakeDomain({{"cubic_box_length", "16"}}),
    {{"particle_type", "../forcefield/data.lj"}}));
  mc.add(Potential(MakeLennardJones()));
  mc.add(Potential(MakeModelHardShape(MakeSlabSine(
    MakeFormulaSineWave({{"amplitude", "2"}, {"width", "8"}}),
    { {"dimension", "0"}, {"wave_dimension", "1"}, {"average_bound0", "-5"},
      {"average_bound1", "5"}}))));
  mc.add(MakeMetropolis({{"beta", "0.1"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "2."}}));
  SeekNumParticles(500).with_trial_add().run(&mc);
  const int steps_per = 1e2;
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)},
                          {"file_name", "tmp/sine"}}));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(1e3);
}

}  // namespace feasst

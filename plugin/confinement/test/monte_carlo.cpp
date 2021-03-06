#include "utils/test/utils.h"
#include "utils/include/checkpoint.h"
#include "utils/include/progress_report.h"
#include "math/include/table.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "configuration/include/file_xyz.h"
#include "system/include/lennard_jones.h"
#include "system/include/hard_sphere.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/visit_model_bond.h"
#include "system/include/long_range_corrections.h"
#include "shape/include/sphere.h"
#include "shape/include/slab.h"
#include "shape/include/slab_sine.h"
#include "shape/include/shape_union.h"
#include "shape/include/shape_intersect.h"
#include "shape/include/cylinder.h"
#include "shape/include/half_space_sine.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/check_energy_and_tune.h"
#include "confinement/include/model_hard_shape.h"
#include "confinement/include/model_table_cartesian.h"
#include "confinement/include/always_reject.h"
#include "confinement/include/henry_coefficient.h"
#include "confinement/include/trial_anywhere.h"
#include "ewald/include/ewald.h"
#include "ewald/include/charge_screened.h"
#include "ewald/include/charge_screened_intra.h"
#include "ewald/include/charge_self.h"

namespace feasst {

TEST(MonteCarlo, ShapeUnion) {
  MonteCarlo mc;
  mc.add(Configuration(MakeDomain({{"cubic_box_length", "8"}}),
    {{"particle_type", "../forcefield/data.lj"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeModelHardShape(MakeShapeUnion(
    MakeSphere(
      {{"radius", "2"}},
      Position({{"x", "0"}, {"y", "0"}, {"z", "0"}})),
    MakeSlab({{"dimension", "2"}, {"bound0", "-1"}, {"bound1", "1"}})))));
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
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
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeModelHardShape(porous_network())));
  mc.set(MakeThermoParams({{"beta", "1.5"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "2."}}));
  SeekNumParticles(500).with_trial_add().run(&mc);
  const int steps_per = 1e3;
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)},
                          {"file_name", "tmp/confine"}}));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(1e4);
}

TEST(MonteCarlo, ShapeTable_LONG) {
  MonteCarlo mc;
  auto domain = MakeDomain({{"cubic_box_length", "20"}});
  mc.add(Configuration(domain, {{"particle_type", "../forcefield/data.lj"}}));
  mc.add(MakePotential(MakeLennardJones()));
  auto pore = porous_network();
  mc.add(MakePotential(MakeModelHardShape(pore)));
  const bool read_table = false;
  //const bool read_table = true;
  std::shared_ptr<ModelTableCart3DIntegr> hamaker;
  if (read_table) {
    auto table = MakeTable3D();
    MakeCheckpoint({{"file_name", "tmp/table"}})->read(table.get());
    hamaker = MakeModelTableCart3DIntegr(table);
  } else {
    hamaker = MakeModelTableCart3DIntegr(MakeTable3D({
      {"num0", "101"},
      {"num1", "101"},
      {"num2", "101"},
      {"default_value", "0."}}));
    #ifdef _OPENMP
    hamaker->compute_table_omp(
    #else // _OPENMP
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
  mc.add(MakePotential(hamaker));

//  EXPECT_NEAR(table2->linear_interpolation(0, 0, 0), 0., NEAR_ZERO);
//  mc.add(MakePotential(MakeModelTableCart3DIntegr(table2)));
//  //mc.add(MakePotential(MakeModelTableCart3DIntegr(table)));
//  mc.set(MakeMetropolis({{"beta", "1.5"}, {"chemical_potential", "1."}}));
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
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeModelHardShape(MakeSlabSine(
    MakeFormulaSineWave({{"amplitude", "2"}, {"width", "8"}}),
    { {"dimension", "0"}, {"wave_dimension", "1"}, {"average_bound0", "-5"},
      {"average_bound1", "5"}}))));
  mc.set(MakeThermoParams({{"beta", "0.1"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "2."}}));
  SeekNumParticles(500).with_trial_add().run(&mc);
  const int steps_per = 1e2;
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)},
                          {"file_name", "tmp/sine"}}));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(1e3);
}

System slab(const int num0 = 0, const int num1 = 0, const int num2 = 0) {
  System system;
  Configuration config(
    MakeDomain({{"cubic_box_length", "5"}}), {
      {"particle_type0", "../forcefield/data.dimer"},
      {"particle_type1", install_dir() + "/plugin/confinement/forcefield/data.slab5x5"},
      {"particle_type2", "../forcefield/data.lj"}});
  EXPECT_EQ(3, config.num_site_types());
  for (int site_type = 0; site_type < config.num_site_types(); ++site_type) {
    config.set_model_param("cutoff", site_type, 2.5);
  }
  for (int i = 0; i < num0; ++i) config.add_particle_of_type(0);
  for (int i = 0; i < num1; ++i) config.add_particle_of_type(1);
  for (int i = 0; i < num2; ++i) config.add_particle_of_type(2);
  system.add(config);
  system.add(MakePotential(MakeLennardJones()));
//  system.add(MakePotential(MakeModelHardShape(MakeSlab({
//    {"dimension", "2"},
//    {"bound0", "2"},
//    {"bound1", "-2"}}))));
  return system;
}

Accumulator henry(System system) {
  MonteCarlo mc;
  //WARN("remove seed");
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.set(system);
  mc.set(MakeThermoParams({{"beta", "1.0"}, {"chemical_potential0", "1"}}));
  mc.set(MakeAlwaysReject());
  mc.add(MakeTrialAdd({{"particle_type", "0"}, {"new_only", "true"}}));
  //mc.add(MakeLogAndMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/henry"}}));
  const int henry_index = mc.num_analyzers();
  mc.add(MakeHenryCoefficient());
  mc.attempt(1e6);
//  const double en_bare = mc.criteria().current_energy();
//  INFO("en_bare " << en_bare);
//  INFO("num " << mc.configuration().num_particles());
//  INFO("parts: " << mc.configuration().selection_of_all().str());
//  for (int i = 0; i < 1e6; ++i) {
//    mc.attempt(1);
//    INFO("en: " << mc.trial(0).accept().energy_new());
//  }
  return mc.analyze(henry_index).accumulator();
}

TEST(ModelTableCart3DIntegr, atomistic_slab_henry_LONG) {
  const Accumulator h = henry(slab(0, 1));
  EXPECT_NEAR(h.average(), 55.5, 4*h.stdev_of_av());
}

TEST(ModelTableCart3DIntegr, table_slab_henry_LONG) {
  System table_system = slab(0, 1, 1);
  table_system.energy();
  auto table = MakeTable3D({
    {"num0", "51"},
    {"num1", "51"},
    {"num2", "51"},
    {"default_value", "0."}});
  auto model = MakeModelTableCart3DIntegr(table);
  const Particle moving_particle = table_system.configuration().particle(1);
  EXPECT_EQ(moving_particle.type(), 2);
  Select select(1, moving_particle);
  #ifdef _OPENMP
    model->compute_table_omp(&table_system, &select);
  #else // _OPENMP
    model->compute_table(&table_system, &select);
  #endif // _OPENMP
  System system = slab(0);
  system.add(MakePotential(model));  // use table instead of explicit wall
  const Accumulator h = henry(system);
  EXPECT_NEAR(h.average(), 52.5, 3*h.stdev_of_av());
}

//TEST(ModelTableCart3DIntegr, atomistic_slab_LONG) {
//  System system = slab(0, 1, 1);
//  system.add(MakePotential(MakeLennardJones()));
//  system.add(MakePotential(MakeModelHardShape(MakeSlab({
//    {"dimension", "2"},
//    {"bound0", "2"},
//    {"bound1", "-2"}}))));
//  auto model_table = MakeModelTableCart3DIntegr(MakeTable3D({
//    {"num0", "51"},
//    {"num1", "51"},
//    {"num2", "51"},
//    {"default_value", "0."}}));
//  Select select(0, system.configuration().particle(0));
//  #ifdef _OPENMP
//    model_table->compute_table_omp(&system, &select);
//  #else // _OPENMP
//    model_table->compute_table(&system, &select);
//  #endif // _OPENMP
//  MakeCheckpoint({{"file_name", "tmp/table2"}})->write(model_table->table());
//  system.add(MakePotential(model_table));
//  // With tabular potential, no longer need data.slab
//  system.get_configuration()->remove_particles(system.configuration().selection_of_all());
//
//  MonteCarlo mc;
//  mc.set(system);
//  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
//  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
//  mc.add(MakeLogAndMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/slabtab"}}));
//  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e4)}, {"tolerance", str(1e-9)}}));
//  SeekNumParticles(15).with_trial_add().run(&mc);
//  mc.attempt(1e6);
//}

TEST(Ewald, henry_coefficient_LONG) {
  System system;
  Configuration config(
    MakeDomain({{"cubic_box_length", "20"}}), {
      {"particle_type0", "../forcefield/data.spce"},
      //{"particle_type0", "../plugin/ewald/forcefield/data.rpm_plus"},
      {"particle_type1", install_dir() + "/plugin/confinement/forcefield/data.slab20x20"}});
  for (int site_type = 0; site_type < config.num_site_types(); ++site_type) {
    config.set_model_param("cutoff", site_type, 2.5);
  }
  config.add_particle_of_type(1);
  FileXYZ().write_for_vmd("tmp.xyz", config);
  system.add(config);
  system.add(MakePotential(
    MakeEwald({{"kmax_squared", "38"},
               {"alpha", str(5.6/system.configuration().domain().min_side_length())}})));
  system.add(MakePotential(MakeModelTwoBodyFactory({MakeLennardJones(),
                                                    MakeChargeScreened()})));
  system.add(MakePotential(MakeChargeScreenedIntra(), MakeVisitModelBond()));
  system.add(MakePotential(MakeChargeSelf()));
  system.add(MakePotential(MakeLongRangeCorrections()));
  DEBUG(system.energy());
  Accumulator h = henry(system);
  DEBUG(h.str());
  EXPECT_NEAR(h.average(), 0.75, 0.03);
}

TEST(HardShape, henry_LONG) {
  for (const double length : {10, 20}) {
    System system;
    system.add(Configuration(MakeDomain({{"cubic_box_length", str(length)}}),
      {{"particle_type0", "../forcefield/data.hard_sphere"}}));
    system.add(MakePotential(MakeModelHardShape(MakeSlab({
      {"dimension", "2"},
      {"bound0", "3"},
      {"bound1", "-3"}}))));
    Accumulator h = henry(system);
    INFO(h.str());
    EXPECT_NEAR(h.average(), 5/system.configuration().domain().min_side_length(), 5*h.block_stdev());
  }
}

}  // namespace feasst

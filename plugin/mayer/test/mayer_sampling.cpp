#include <cmath>
#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/utils.h"
#include "system/include/hard_sphere.h"
#include "system/include/lennard_jones.h"
#include "system/include/model_two_body_factory.h"
#include "models/include/square_well.h"
#include "steppers/include/log.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/movie.h"
#include "steppers/include/tune.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/run.h"
#include "mayer/include/mayer_sampling.h"
#include "ewald/include/coulomb.h"
#include "ewald/include/utils.h"
#include "models/include/lennard_jones_force_shift.h"

namespace feasst {

TEST(MayerSampling, ljb2) {
  System system = two_particle_system();
  system.add_to_reference(MakePotential(MakeHardSphere()));
  auto translate = MakeTrialTranslate({{"new_only", "true"}, {"reference_index", "0"}, {"weight", "0.75"}});
  //auto translate = MakeTrialTranslate({{"tunable_param", "0.5"}});
  /// HWH notes: does this need a max?
  const int nTrialsEq = 1e4, nTrials = 2e4;
  //const int nTrialsEq = 1e6, nTrials = 1e6;
  Configuration * config = system.get_configuration();
  config->set_model_param("cutoff", 0, NEAR_INFINITY);
  EXPECT_EQ(config->model_params().cutoff().value(0), NEAR_INFINITY);
  const double boxl = 2*(config->model_params().cutoff().value(0));
  config->set_side_lengths(Position().set_vector({boxl, boxl, boxl}));
  std::cout << "boxl " << boxl << std::endl;
  system.set(MakeThermoParams({{"beta", "1."},
    {"chemical_potential", "-2.775"}}));
  MayerSampling criteria;
  criteria.set_current_energy(system.energy());
  RandomMT19937 random;
  for (int iTrial = 0; iTrial < nTrialsEq + nTrials; ++iTrial) {
    translate->attempt(&criteria, &system, &random);
  }
  const double b2 = 2./3.*PI*criteria.second_virial_ratio();
  std::cout << "b2 " << b2 << std::endl;
  EXPECT_NEAR(-5.3, b2, 15);
  EXPECT_GT(std::abs(2.0944-b2), 0.0001); // HS value

  std::shared_ptr<Criteria> crit2 = test_serialize<MayerSampling, Criteria>(criteria);
  EXPECT_EQ(system.thermo_params().beta(), 1.);
}

MayerSampling ljb2(const int trials) {
  MonteCarlo mc;
  { // initialize system
    System lj = lennard_jones({{"cubic_box_length", "1000"}});
    Configuration * config = lj.get_configuration();
    config->set_model_param("cutoff", 0, config->domain().side_length(0)/2.);
    config->add_particle_of_type(0);
    config->add_particle_of_type(0);
    mc.set(lj);
  }
  mc.add_to_reference(MakePotential(MakeHardSphere()));
  mc.set(MakeThermoParams({{"beta", "1."}, {"chemical_potential", "-2.775"}}));
  mc.set(MakeMayerSampling());
  mc.add(MakeTrialTranslate({{"new_only", "true"}, {"reference_index", "0"}, {"weight", "0.75"}}));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(trials);
  std::stringstream ss;
  mc2.criteria().serialize(ss);
  MayerSampling mayer(ss);
  return mayer;
}

TEST(MonteCarlo, ljb2) {
  MayerSampling mayer = ljb2(1e4);
  const double b2 = 2./3.*PI*mayer.second_virial_ratio();
  INFO("b2 " << b2);
  EXPECT_NEAR(-5.3, b2, 20);
  EXPECT_GT(std::abs(2.0944 - b2), 0.0001); // HS value
}

TEST(MonteCarlo, ljb2_LONG) {
  MayerSampling mayer = ljb2(1e7);
  const double b2 = 2./3.*PI*mayer.second_virial_ratio();
  INFO("b2 " << b2);
  EXPECT_NEAR(-5.3, b2, 0.3);
  EXPECT_GT(std::abs(2.0944 - b2), 0.0001); // HS value
}

// Check SPCE

TEST(MayerSampling, cg4_rigid_LONG) {
  MonteCarlo mc;
  auto config = MakeConfiguration({{"cubic_box_length", "1000"},
    {"particle_type0", install_dir() + "/plugin/chain/forcefield/cg4_mab.fstprt"}});
  config->add_particle_type(install_dir() + "/plugin/chain/forcefield/cg4_mab.fstprt", "2");
  config->add_particle_of_type(0);
  config->add_particle_of_type(1);
  mc.add(config);
  EXPECT_EQ(2, mc.configuration().num_particles());
  EXPECT_EQ(1, mc.configuration().num_particles_of_type(0));
  mc.add(MakePotential(MakeSquareWell()));
  mc.add_to_reference(MakePotential(MakeHardSphere()));
  const double temperature = 0.7092;
  mc.set(MakeThermoParams({{"beta", str(1./temperature)}}));
  mc.set(MakeMayerSampling());
  mc.add(MakeTrialTranslate({{"new_only", "true"}, {"reference_index", "0"},
    {"tunable_param", "1"}, {"particle_type", "1"}}));
  mc.add(MakeTrialRotate({{"new_only", "true"}, {"reference_index", "0"},
    {"tunable_param", "40"}}));
  std::string steps_per = "1e4";
  mc.add(MakeLog({{"steps_per", steps_per}, {"file_name", "tmp/cg4.txt"}}));
  mc.add(MakeMovie({{"steps_per", steps_per}, {"file_name", "tmp/cg4.xyz"}}));
  mc.add(MakeTune({{"steps_per", steps_per}}));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(1e6);
  std::stringstream ss;
  mc2.criteria().serialize(ss);
  MayerSampling mayer(ss);
  EXPECT_NEAR(0.46, mayer.second_virial_ratio(), 0.15);
}

//// HWH: there is something wrong with mayer sampling + Coulomb
//// Table 2 of https://pubs.acs.org/doi/pdf/10.1021/jp0710685
//// Fig 1 of https://doi.org/10.1063/1.5016165
//TEST(MayerSampling, SPCE_LONG) {
//  MonteCarlo mc;
//  { auto config = MakeConfiguration({{"cubic_box_length", str(NEAR_INFINITY)}});
//    config->add_particle_type(install_dir() + "/forcefield/spce.fstprt");
//    config->add_particle_type(install_dir() + "/forcefield/spce.fstprt", "2");
//    for (int stype = 0; stype < config->num_site_types(); ++stype) {
//      config->set_model_param("cutoff", stype, 100.);
//      // HWH Why is there dependence on the cutoff? If large, it drifts too far away
//      //config->set_model_param("cutoff", stype, config->domain().side_length(0)/2.);
//    }
//    config->add_particle_of_type(0);
////    config->add_particle_of_type(0);
//    config->add_particle_of_type(1);
//    mc.add(config);
//  }
//  mc.add(MakePotential(MakeModelTwoBodyFactory({MakeLennardJones(), MakeCoulomb()})));
//  //mc.add_to_reference(MakePotential(MakeLennardJones()));
//  auto ref = MakePotential(MakeHardSphere());
//  ref->set_model_params(mc.configuration());
//  ref->set_model_param("sigma", 0, 30);
//  ref->set_model_param("sigma", 2, 30);
//  mc.add_to_reference(ref);
//  //const double temperature = 373; // kelvin
//  //const double temperature = 300; // kelvin
//  //const double temperature = 400; // kelvin
//  const double temperature = 500; // kelvin
//  //const double temperature = 1e3; // kelvin
//  //const double temperature = 773; // kelvin
//  mc.set(MakeThermoParams({{"beta", str(1./kelvin2kJpermol(temperature))}}));
//  auto mayer = MakeMayerSampling();
//  mc.set(mayer);
//  //mc.add(MakeTrialTranslate({{"new_only", "true"}, {"reference_index", "0"},
//  //  {"tunable_param", "1"}}));
//  mc.add(MakeTrialTranslate({{"new_only", "true"}, {"reference_index", "0"},
//    {"tunable_param", "1"}, {"particle_type", "1"}}));
//  mc.add(MakeTrialRotate({{"new_only", "true"}, {"reference_index", "0"},
//    {"tunable_param", "40"}}));
//  std::string steps_per = "1e5";
//  mc.add(MakeLog({{"steps_per", steps_per}, {"file_name", "tmp/spce.txt"}}));
//  mc.add(MakeMovie({{"steps_per", steps_per}, {"file_name", "tmp/spce.xyz"}}));
//  mc.add(MakeTune({{"steps_per", steps_per}}));
//  mc.attempt(1e6);
//  mc.perform(MakeRemoveModify({{"name", "Tune"}}));
//  mc.attempt(1e7);
//  mayer = MakeMayerSampling();
//  mc.set(mayer);
//  mc.attempt(1e7);
//  INFO("mayer ref " << mayer->mayer_ref().str());
//  double b2hs = 2./3.*PI*std::pow(mc.system().reference(0, 0).model_params().sigma().value(0), 3); // A^3
//  INFO(mc.system().reference(0, 0).model_params().sigma().value(0));
//  //double b2hs = 2./3.*PI*std::pow(mc.configuration().model_params().sigma().value(0), 3); // A^3
//  INFO("b2hs(A^3) " << b2hs);
//  b2hs *= 1e-30*1e3*mc.configuration().physical_constants().avogadro_constant(); // L/mol
//  INFO("b2hs(L/mol) " << b2hs);
//  INFO("b2spce/b2hs " << mayer->second_virial_ratio());
//  INFO("b2spce(L/mol) " << b2hs*mayer->second_virial_ratio());
//  EXPECT_NEAR(-0.4082, b2hs*mayer->second_virial_ratio(), 0.006);
//  //EXPECT_NEAR(-0.08596, b2hs*mayer->second_virial_ratio(), 0.006);
//}

// https://dx.doi.org/10.1063/1.4918557
TEST(MayerSampling, trimer_LONG) {
  MonteCarlo mc;
  { auto config = MakeConfiguration({{"cubic_box_length", str(NEAR_INFINITY)}});
    config->add_particle_type(install_dir() + "/forcefield/trimer_0.4L.fstprt");
    config->add_particle_type(install_dir() + "/forcefield/trimer_0.4L.fstprt", "2");
    config->add_particle_of_type(0);
    config->add_particle_of_type(1);
    const double rwca = std::pow(2, 1./6.);
    config->set_model_param("cutoff", 0, 1, rwca);
    config->set_model_param("cutoff", 0, 3, rwca);
    config->set_model_param("cutoff", 1, 2, rwca);
    config->set_model_param("cutoff", 2, 3, rwca);
    EXPECT_NEAR(config->model_params().mixed_cutoff()[0][0], 3, 1e-14);
    EXPECT_NEAR(config->model_params().mixed_cutoff()[0][1], rwca, 1e-14);
    EXPECT_NEAR(config->model_params().mixed_cutoff()[0][2], 3, 1e-14);
    EXPECT_NEAR(config->model_params().mixed_cutoff()[0][3], rwca, 1e-14);
    EXPECT_NEAR(config->model_params().mixed_cutoff()[1][1], rwca, 1e-14);
    EXPECT_NEAR(config->model_params().mixed_cutoff()[1][2], rwca, 1e-14);
    EXPECT_NEAR(config->model_params().mixed_cutoff()[1][3], rwca, 1e-14);
    EXPECT_NEAR(config->model_params().mixed_cutoff()[2][2], 3, 1e-14);
    EXPECT_NEAR(config->model_params().mixed_cutoff()[2][3], rwca, 1e-14);
    EXPECT_NEAR(config->model_params().mixed_cutoff()[3][3], rwca, 1e-14);
    mc.add(config);
  }
  mc.add(MakePotential(MakeLennardJonesForceShift()));
  auto ref = MakePotential(MakeHardSphere());
  auto params = ref->model_params(mc.system().configuration());
  for (int i = 0; i < mc.system().configuration().num_site_types(); ++i) {
    for (int j = 0; j < mc.system().configuration().num_site_types(); ++j) {
      params.set("sigma", i, j, 0.);
    }
  }
  params.set("sigma", 0, 0, 1.);
  params.set("sigma", 0, 2, 1.);
  params.set("sigma", 2, 2, 1.);
  ref->set(params);
  mc.add_to_reference(ref);
  mc.set(MakeThermoParams({{"beta", str(1./0.815)}}));
  auto mayer = MakeMayerSampling();
  mc.set(mayer);
  mc.add(MakeTrialTranslate({{"new_only", "true"}, {"reference_index", "0"},
    {"tunable_param", "1"}, {"particle_type", "1"}}));
  mc.add(MakeTrialRotate({{"new_only", "true"}, {"reference_index", "0"},
    {"tunable_param", "40"}}));
  const std::string steps_per = "1e4";
  mc.add(MakeLogAndMovie({{"steps_per", steps_per}, {"file_name", "tmp/trib"}}));
  mc.attempt(1e6);
  double b2hs = 2./3.*PI*std::pow(mc.configuration().model_params().sigma().value(0), 3); // A^3
  INFO(mayer->second_virial_ratio());
  INFO(b2hs*mayer->second_virial_ratio());
  INFO("mayer: " << mayer->mayer().str());
  INFO("mayer_ref: " << mayer->mayer_ref().str());
  EXPECT_NEAR(0, mayer->mayer().average(), 4*mayer->mayer().block_stdev());
}

}  // namespace feasst

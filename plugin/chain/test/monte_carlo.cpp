#include <memory>
#include "utils/test/utils.h"
#include "utils/include/utils_io.h"
#include "math/include/accumulator.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "monte_carlo/include/metropolis.h"
#include "system/test/system_test.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/check_physicality.h"
#include "chain/include/analyze_rigid_bonds.h"
#include "chain/include/trial_grow.h"
#include "chain/include/trial_pivot.h"
#include "chain/include/trial_crankshaft.h"
#include "chain/include/trial_reptate.h"
#include "chain/include/trial_swap_sites.h"
#include "chain/include/utils_chain.h"
#include "chain/include/recenter_particles.h"
#include "ewald/include/ewald.h"
#include "ewald/test/system_example.h"

namespace feasst {

Configuration config() {
  Configuration config({
    {"cubic_box_length", "12"},
    {"particle_type0", "../forcefield/data.chain10"},
    //{"particle_type0", "../plugin/chain/forcefield/data.chain50"},
    {"init_cells", "1."},
  });
  config.add_particle_of_type(0);
  return config;
}

Potential lj_dual_cut(const Configuration config) {
  Potential lj_dual_cut(MakeLennardJones(), MakeVisitModelCell());
  lj_dual_cut.set_model_params(config);
  lj_dual_cut.set_model_param("cutoff", 0, 1);
  return lj_dual_cut;
}

Potential lj_intra_dual_cut(const Configuration config) {
  Potential lj_intra_dual_cut(MakeLennardJones(),
                              MakeVisitModelIntra({{"cutoff", "1"}}));
  lj_intra_dual_cut.set_model_params(config);
  lj_intra_dual_cut.set_model_param("cutoff", 0, 1);
  return lj_intra_dual_cut;
}

System chain_system() {
  System system;
  system.add(config());
  system.add_to_unoptimized(Potential(MakeLennardJones()));
  system.add_to_reference(lj_dual_cut(system.configuration()));
  system.add_to_unoptimized(Potential(MakeLennardJones(),
                                      MakeVisitModelIntra({{"cutoff", "1"}})));
  system.add_to_reference(lj_intra_dual_cut(system.configuration()));
  system.add_to_unoptimized(Potential(MakeLongRangeCorrections()));
  return system;
}

TEST(MonteCarlo, chain) {
  MonteCarlo mc;
  mc.set(chain_system());
  mc.set(MakeMetropolis({{"beta", "1"}, {"chemical_potential", "1."}}));
  mc.seek_num_particles(1);
  mc.add(MakeTrialTranslate({
    {"weight", "1."},
    {"tunable_param", "1."},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialRotate({
    {"weight", "1."},
    {"tunable_param", "20."},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialPivot({
    {"weight", "1."},
    {"tunable_param", "20."},
    {"max_length", "30"},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialReptate({
    {"weight", "1."},
    {"max_length", "1"},
//    {"reference_index", "0"},
//    {"num_steps", "2"},
  }));
  mc.add(MakeTrialCrankshaft({
    {"weight", "1."},
    {"tunable_param", "25."},
    {"max_length", "5."},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialGrowLinear(
    MakeTrialComputeMove(),
    {
//      {"weight", "0.1"},
      {"particle_type", "0"},
      {"num_steps", "3"},
      {"reference_index", "0"},
    }
  ));
  const int steps_per = 1e2;
  mc.add(MakeLog({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/chainlog.txt"},
  }));
  mc.add(MakeMovie({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/chain10movie.xyz"},
  }));
  mc.add(MakeCheckEnergy({
    {"steps_per", str(steps_per)},
    {"tolerance", "1e-10"},
  }));
  mc.add(MakeTuner({{"steps_per", str(steps_per)}}));
  mc.add(MakeAnalyzeRigidBonds({{"steps_per", str(steps_per)}}));
  mc.attempt(3e2);

  MonteCarlo mc2 = test_serialize(mc);
  EXPECT_EQ(mc2.analyzers().size(), 3);
}

TEST(MonteCarlo, deprotonation) {
  MonteCarlo monte_carlo;
  monte_carlo.set(MakeRandomMT19937({{"seed", "time"}}));
  //monte_carlo.set(MakeRandomMT19937({{"seed", "1582203788"}}));
  { Configuration config({
      {"cubic_box_length", "20"},
      {"particle_type0", "../plugin/chain/forcefield/data.chain10titratable"},
      {"particle_type1", "../forcefield/data.lj"}
    });
    config.add_particle_of_type(0);
    monte_carlo.add(config);
  }
  monte_carlo.add(Potential(MakeLennardJones()));
  monte_carlo.add(Potential(MakeLennardJones(), MakeVisitModelIntra({{"cutoff", "1"}})));
  monte_carlo.set(MakeMetropolis({
    {"beta", "1"},
    {"chemical_potential0", "1"},
    {"chemical_potential1", "-5"},
    {"pH", "0.001"} // HWH something appears to be wrong with pH definition
  }));

  // trials for chains
  monte_carlo.add(MakeTrialTranslate({{"weight", "1."}, {"particle_type", "0"}, {"tunable_param", "2"}}));
  monte_carlo.add(MakeTrialRotate({{"weight", "1."}, {"particle_type", "0"}, {"tunable_param", "50."}}));
  monte_carlo.add(MakeTrialPivot({{"weight", "1."}, {"particle_type", "0"}, {"tunable_param", "50."}}));
  monte_carlo.add(MakeTrialCrankshaft({{"weight", "1."}, {"particle_type", "0"}, {"tunable_param", "50."}}));
  monte_carlo.add(MakeTrialReptate({{"weight", "1."}, {"particle_type", "0"}, {"max_length", "1"}}));
  monte_carlo.add(MakeTrialSwapSites({
      {"weight", "1."},
      {"particle_type", "0"},
      {"site_type1", "0"},
      {"site_type2", "1"}
  }));
  monte_carlo.add(MakeTrialGrowLinear(MakeTrialComputeMove(), {
      {"weight", str(1./4./10.)},
      {"particle_type", "0"},
      {"num_steps", "4"}
  }));
  add_deprotonation(&monte_carlo, {
      {"reactant_type", "0"},
      {"reactant_site_type", "0"},
      {"new_site_type", "1"},
      {"add_type", "1"},
  });

  // trials for counterions
  monte_carlo.add(MakeTrialTranslate({{"weight", "1."}, {"particle_type", "1"}, {"tunable_param", "2"}}));
  add_trial_transfer(&monte_carlo, {{"particle_type", "1"}});

  // analysis
  const int steps_per = 1e3;
  monte_carlo.add(MakeCheckEnergy({{"steps_per", str(steps_per)},
                                   {"tolerance", str(1e-10)}}));
  monte_carlo.add(MakeLog({{"steps_per", str(steps_per)}, {"file_name", "tmp/titra.txt"},
                           {"clear_file", "true"}}));
  monte_carlo.add(MakeMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/titra.xyz"},
                             {"clear_file", "true"}}));
  monte_carlo.add(MakeTuner({{"steps_per", str(steps_per)}}));
  monte_carlo.add(MakeAnalyzeRigidBonds({{"steps_per", str(steps_per)}}));
  EXPECT_EQ(1, monte_carlo.configuration().particle_type_to_group(0));
  monte_carlo.add(MakeRecenterParticles({{"steps_per", str(steps_per)},
    {"group_index", str(monte_carlo.configuration().particle_type_to_group(0))}}));
  monte_carlo.add(MakeCheckPhysicality({{"steps_per", str(steps_per)}}));

  // perform simulation
  monte_carlo.attempt(1e4);
}

TEST(MonteCarlo, deprotonationWithEwald) {
  MonteCarlo monte_carlo;
  monte_carlo.set(MakeRandomMT19937({{"seed", "time"}}));
  //monte_carlo.set(MakeRandomMT19937({{"seed", "1582203788"}}));
  { System system;
    { Configuration config({
        {"cubic_box_length", "20"},
        {"particle_type0", "../plugin/chain/forcefield/data.chain10titratable"},
        {"particle_type1", "../plugin/ewald/forcefield/data.rpm_plus"},
        {"particle_type2", "../plugin/ewald/forcefield/data.rpm_minus"},
      });
      config.add_particle_of_type(0);
      system.add(config);
    }
    system.add(Potential(
      MakeEwald({{"kmax_squared", "38"},
                 {"alpha", str(5.6/system.configuration().domain().min_side_length())}}),
      {{"prevent_cache", "true"}}));
    system.add(Potential(MakeModelTwoBodyFactory({MakeLennardJones(),
                                               MakeChargeScreened()})));
    system.add(Potential(MakeChargeScreenedIntra(),
                         MakeVisitModelIntra({{"cutoff", "0"}})));
    system.add(Potential(MakeChargeSelf()));
    system.add(Potential(MakeLongRangeCorrections()));
    //system.add(Potential(MakeLennardJones(), MakeVisitModelIntra({{"cutoff", "1"}})));
    monte_carlo.set(system);
  }
  monte_carlo.set(MakeMetropolis({
    {"beta", "1"},
    {"chemical_potential0", "1"},
    {"chemical_potential1", "-5"},
    {"pH", "0.001"} // HWH something appears to be wrong with pH definition
  }));

  // trials for chains
  monte_carlo.add(MakeTrialTranslate({{"weight", "1."}, {"particle_type", "0"}, {"tunable_param", "2"}}));
  monte_carlo.add(MakeTrialRotate({{"weight", "1."}, {"particle_type", "0"}, {"tunable_param", "50."}}));
  monte_carlo.add(MakeTrialPivot({{"weight", "1."}, {"particle_type", "0"}, {"tunable_param", "50."}}));
//  monte_carlo.add(MakeTrialCrankshaft({{"weight", "1."}, {"particle_type", "0"}, {"tunable_param", "50."}}));
//  monte_carlo.add(MakeTrialReptate({{"weight", "1."}, {"particle_type", "0"}, {"max_length", "1"}}));
//  monte_carlo.add(MakeTrialSwapSites({
//      {"weight", "1."},
//      {"particle_type", "0"},
//      {"site_type1", "0"},
//      {"site_type2", "1"}
//  }));
//  monte_carlo.add(MakeTrialGrowLinear(MakeTrialComputeMove(), {
//      {"weight", str(1./4./10.)},
//      {"particle_type", "0"},
//      {"num_steps", "4"}
//  }));
//  add_deprotonation(&monte_carlo, {
//      {"reactant_type", "0"},
//      {"reactant_site_type", "0"},
//      {"new_site_type", "1"},
//      {"add_type", "1"},
//  });

  // trials for counterions
//  monte_carlo.add(MakeTrialTranslate({{"weight", "1."}, {"particle_type", "1"}, {"tunable_param", "2"}}));
//  add_trial_transfer(&monte_carlo, {{"particle_type", "1"}});

  // analysis
  const int steps_per = 1e3;
  monte_carlo.add(MakeCheckEnergy({{"steps_per", str(steps_per)},
                                   {"tolerance", str(1e-10)}}));
  monte_carlo.add(MakeLog({{"steps_per", str(steps_per)}, {"file_name", "titra.txt"},
                           {"clear_file", "true"}}));
  monte_carlo.add(MakeMovie({{"steps_per", str(steps_per)}, {"file_name", "titra.xyz"},
                             {"clear_file", "true"}}));
  monte_carlo.add(MakeTuner({{"steps_per", str(steps_per)}}));
  monte_carlo.add(MakeAnalyzeRigidBonds({{"steps_per", str(steps_per)}}));
  EXPECT_EQ(1, monte_carlo.configuration().particle_type_to_group(0));
  monte_carlo.add(MakeRecenterParticles({{"steps_per", str(steps_per)},
    {"group_index", str(monte_carlo.configuration().particle_type_to_group(0))}}));
  monte_carlo.add(MakeCheckPhysicality({{"steps_per", str(steps_per)}}));

  // perform simulation
  monte_carlo.attempt(1e3);
}

}  // namespace feasst

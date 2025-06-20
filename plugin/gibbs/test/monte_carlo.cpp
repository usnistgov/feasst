#include <memory>
#include <cmath>
#include "utils/test/utils.h"
#include "utils/include/io.h"
#include "utils/include/checkpoint.h"
#include "utils/include/progress_report.h"
#include "math/include/accumulator.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/potential.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"
#include "system/include/dont_visit_model.h"
#include "system/include/ideal_gas.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/constrain_num_particles.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_volume.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/run.h"
#include "steppers/include/num_particles.h"
#include "steppers/include/movie.h"
#include "steppers/include/tune.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/profile_trials.h"
#include "gibbs/include/trial_gibbs_particle_transfer.h"
#include "gibbs/include/trial_gibbs_volume_transfer.h"

namespace feasst {

TEST(MonteCarlo, gibbs_ensemble) {
  auto trial = MakeTrialGibbsParticleTransfer();
  auto mc = MakeMonteCarlo({{
    //{"RandomMT19937", {{"seed", "12345"}}},
    {"Configuration", {{"name", "first"}, {"xyz_file", "../plugin/configuration/test/data/lj_sample_config_periodic4.xyz"},
      {"particle_type", "0:../particle/lj.txt"}, {"group0", "first"}, {"first_particle_index", "0"}, {"cutoff", "2"}}},
    {"Configuration", {{"name", "second"}, {"xyz_file", "../plugin/configuration/test/data/lj_sample_config_periodic4.xyz"},
      {"particle_type", "0:../particle/lj.txt"}, {"group0", "first"}, {"first_particle_index", "0"}, {"cutoff", "2"}}},
//    {"CopyFollowingLines", {{"for_num_configurations", "2"}}},
    {"Potential", {{"Model", "LennardJones"}, {"config", "first"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}, {"config", "first"}}},
    {"RefPotential", {{"VisitModel", "DontVisitModel"}, {"config", "first"}, {"ref", "first"}}},
    {"Potential", {{"Model", "LennardJones"}, {"config", "second"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}, {"config", "second"}}},
    {"RefPotential", {{"VisitModel", "DontVisitModel"}, {"config", "second"}, {"ref", "first"}}},
//    {"EndCopy", {{}}},
    {"ThermoParams", {{"beta", str(1/0.8)}, {"chemical_potential", "1."}, {"pressure", "1"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"config", "first"}}},
    {"TrialTranslate", {{"config", "second"}}},
    {"TrialGrowFile", {{"grow_file", "../plugin/gibbs/test/data/grow.txt"}}},
    {"TrialGibbsParticleTransfer", {{"particle_type", "0"}, {"ref", "first"}, {"configs", "first,second"}}},
    {"TrialGibbsVolumeTransfer", {{"tunable_param", "0.1"}, {"ref", "first"}, {"configs", "first,second"}}},
    {"Tune", {{"trials_per_tune", "4"}}},
    {"Log", {{"trials_per_write", "1e0"}, {"output_file", "tmp/lj.txt"}}},
    {"Movie", {{"trials_per_write", "1e0"}, {"output_file", "tmp/lj0.xyz"}, {"config", "first"}}},
    {"Movie", {{"trials_per_write", "1e0"}, {"output_file", "tmp/lj1.xyz"}, {"config", "second"}}},
    {"CheckEnergy", {{"trials_per_update", "1e0"}, {"tolerance", str(1e-9)}}},
    {"CheckConstantVolume", {{"trials_per_update", "1e0"}}},
    {"GhostTrialVolume", {{"trials_per_update", str(1e0)}, {"trials_per_write", "1"}, {"output_file", "tmp/lj_p.csv"}}},
    //{"Run", {{"num_trials", "1e2"}}},
    //{"Run", {{"num_trials", "1e6"}}},
  }}, true);
  //}});
  auto mc2 = test_serialize_unique(*mc);
  mc2->run_num_trials(1e2);
  EXPECT_EQ(2, mc2->system().num_configurations());
  EXPECT_NE(mc2->system().potential(0, 0).stored_energy(),
            mc2->system().potential(0, 1).stored_energy());
  EXPECT_EQ(mc2->system().configuration(0).num_particles() +
            mc2->system().configuration(1).num_particles(), 60);
  EXPECT_NEAR(mc2->system().configuration(0).domain().volume() +
              mc2->system().configuration(1).domain().volume(),
              2*std::pow(8,3),
              1e-8);
}

TEST(MonteCarlo, gibbs_adjust) {
  auto mc = MakeMonteCarlo({{
    //{"RandomMT19937", {{"seed", "1728330862"}}},
    //{"CopyNextLine", {{}}},
    {"Configuration", {{"xyz_file", "../plugin/configuration/test/data/lj_sample_config_periodic4.xyz"},
      {"particle_type0", "../particle/lj.txt"}, {"cutoff", "2"}}},
    {"Configuration", {{"cubic_side_length", "20"}, {"particle_type0", "../particle/lj.txt"}, {"add_particles_of_type0", "1"}}},
    {"CopyFollowingLines", {{"for_num_configurations", "2"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"RefPotential", {{"VisitModel", "DontVisitModel"}}},
    {"EndCopy", {{}}},
    {"ThermoParams", {{"beta", str(1/0.8)}, {"chemical_potential", "1."}, {"pressure", "1"}}},
    {"Metropolis", {{}}},
    {"CopyNextLine", {{"replace", "configuration_index"}, {"with", "0"}}},
    {"TrialTranslate", {{"configuration_index", "1"}}},
    {"TrialGrowFile", {{"grow_file", "../plugin/gibbs/test/data/grow.txt"}}},
    {"TrialGibbsParticleTransfer", {{"particle_type", "0"}, {"reference_index", "0"}}},
    {"TrialGibbsVolumeTransfer", {{"tunable_param", "0.1"}, {"reference_index", "0"}}},
    {"Tune", {{"trials_per_tune", "4"}}},
    {"Log", {{"trials_per_write", "1e0"}, {"output_file", "tmp/lj.txt"}}},
    {"CopyFollowingLines", {{"for_num_configurations", "2"}, {"replace_with_index", "[config]"}}},
    {"Movie", {{"trials_per_write", "1e0"}, {"output_file", "tmp/ljc[config].xyz"}}},
    {"Density", {{"trials_per_write", "1e0"}, {"output_file", "tmp/ljc[config]dens.xyz"}}},
    {"EndCopy", {{}}},
    {"CheckEnergy", {{"trials_per_update", "1e0"}, {"tolerance", str(1e-9)}}},
    {"GibbsInitialize", {{"trials_per_update", "1e0"}}},
    {"Run", {{"num_trials", "1e2"}}},
  //}});
  }}, true);
  auto mc2 = test_serialize_unique(*mc);
  EXPECT_EQ(2, mc2->system().num_configurations());
}

}  // namespace feasst

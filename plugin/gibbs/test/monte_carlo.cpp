#include <memory>
#include <cmath>
#include "utils/test/utils.h"
#include "utils/include/io.h"
#include "utils/include/checkpoint.h"
#include "utils/include/progress_report.h"
#include "math/include/accumulator.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"
#include "system/include/dont_visit_model.h"
#include "system/include/ideal_gas.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/constrain_num_particles.h"
#include "monte_carlo/include/seek_num_particles.h"
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
    //{"RandomMT19937", {{"seed", "1234"}}},
    //{"RandomMT19937", {{"seed", "1693947707"}}},
    {"Configuration", {{"xyz_file", "../plugin/configuration/test/data/lj_sample_config_periodic4.xyz"},
      {"particle_type0", "../particle/lj.fstprt"}, {"group0", "first"}, {"first_particle_index", "0"}}},
    {"Configuration", {{"xyz_file", "../plugin/configuration/test/data/lj_sample_config_periodic4.xyz"},
      {"particle_type0", "../particle/lj.fstprt"}, {"group0", "first"}, {"first_particle_index", "0"}}},
    {"Potential", {{"Model", "LennardJones"}, {"configuration_index", "0"}}},
    //{"Potential", {{"VisitModel", "LongRangeCorrections"}, {"configuration_index", "0"}}},
    {"Potential", {{"Model", "LennardJones"}, {"configuration_index", "1"}}},
    //{"Potential", {{"VisitModel", "LongRangeCorrections"}, {"configuration_index", "1"}}},
    {"ThermoParams", {{"beta", "100.2"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"configuration_index", "0"}}},
    {"TrialTranslate", {{"configuration_index", "1"}}},
    {"TrialGibbsParticleTransfer", {{"particle_type", "0"}}},
    //{"TrialGibbsVolumeTransfer", {{}}},
    {"TrialGibbsVolumeTransfer", {{"tunable_param", "0.1"}}},
    //{"TrialGibbsParticleTransfer", {{"configuration_index0", "0"}, {"configuration_index1", "1"}}},
    {"Log", {{"trials_per_write", str(1e0)}, {"file_name", "tmp/lj.txt"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"file_name", "tmp/lj0.xyz"}, {"configuration_index", "0"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"file_name", "tmp/lj1.xyz"}, {"configuration_index", "1"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
    {"Run", {{"num_trials", "1e2"}}},
  }});
  EXPECT_EQ(2, mc->system().num_configurations());
  EXPECT_NE(mc->system().potential(0, 0).stored_energy(),
            mc->system().potential(0, 1).stored_energy());
  EXPECT_EQ(mc->system().configuration(0).num_particles() +
            mc->system().configuration(1).num_particles(), 60);
  EXPECT_NEAR(mc->system().configuration(0).domain().volume() +
              mc->system().configuration(1).domain().volume(),
              2*std::pow(8,3),
              1e-8);
}

}  // namespace feasst

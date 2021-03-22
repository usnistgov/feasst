#include "utils/test/utils.h"
#include "math/include/utils_math.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/utils.h"
#include "steppers/include/pair_distribution.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/energy.h"
#include "steppers/include/seek_analyze.h"
#include "steppers/include/seek_modify.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trials.h"

namespace feasst {

TEST(PairDistribution, gr_LONG) {
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  //const int num = 2;
  const int num = 500;
  const double density = 0.776;
  mc.set(lennard_jones({{"cubic_box_length",
    feasst::str(std::pow(500/density, 1./3.))}}));
  mc.set(MakeThermoParams({{"beta", feasst::str(1./0.85)}, {"chemical_potential", "-6."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  //const int num = round(0.776*mc.configuration().domain().volume());
  //const int num = 2;
  SeekNumParticles(num).with_trial_add().run(&mc);
  const std::string steps_per = str(1e5);
  mc.attempt(1e6); // equilibrate
  mc.add(MakeLogAndMovie({{"steps_per", steps_per}, {"file_name", "tmp/lj"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", steps_per}, {"tolerance", str(1e-9)}}));
  mc.add(MakeEnergy({{"steps_per_write", steps_per}, {"file_name", "tmp/ljen.txt"}}));
  mc.add(MakePairDistribution({
    {"steps_per_update", "100"},
    {"steps_per_write", "10000"},
    {"dr", "0.025"},
    {"file_name", "tmp/gr.txt"},
    {"multistate", "true"},
    {"multistate_aggregate", "false"},
  }));
  EXPECT_EQ(mc.configuration().model_params().cutoff().value(0), 3.);
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(1e6);
  MonteCarlo mc3 = test_serialize(mc2);
  auto pd = PairDistribution(mc3.modify(SeekModify().index("PairDistribution", mc3)[0]).modify(0));
  const grtype& gr = pd.radial(mc3.configuration());
  const Accumulator en = Energy(mc3.analyze(SeekAnalyze().index("Energy", mc3)[0])).energy();

  // compare with https://mmlapps.nist.gov/srs/LJ_PURE/md.htm
  EXPECT_EQ(gr[43].first, 1.0875);
  EXPECT_NEAR(gr[43].second[0][0], 2.65, 0.2);
  EXPECT_NEAR(en.average()/num, -5.517, 0.06);
}

}  // namespace feasst

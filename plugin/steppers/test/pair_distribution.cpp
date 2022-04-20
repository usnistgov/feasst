#include "utils/test/utils.h"
#include "math/include/utils_math.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "steppers/include/pair_distribution.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/energy.h"
#include "steppers/include/seek_analyze.h"
#include "steppers/include/seek_modify.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"

namespace feasst {

TEST(PairDistribution, gr_LONG) {
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  //const int num = 2;
  const int num = 500;
  const double density = 0.776;
  mc.add(MakeConfiguration({
    {"cubic_box_length", feasst::str(std::pow(500/density, 1./3.))},
    {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.set(MakeThermoParams({{"beta", feasst::str(1./0.85)}, {"chemical_potential", "-6."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", str(num)}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  const std::string trials_per = str(1e5);
  mc.attempt(1e6); // equilibrate
  mc.add(MakeLogAndMovie({{"trials_per", trials_per}, {"file_name", "tmp/lj"}}));
  mc.add(MakeCheckEnergyAndTune({{"trials_per", trials_per}, {"tolerance", str(1e-9)}}));
  mc.add(MakeEnergy({{"trials_per_write", trials_per}, {"file_name", "tmp/ljen.txt"}}));
  mc.add(MakePairDistribution({
    {"trials_per_update", "100"},
    {"trials_per_write", "10000"},
    {"dr", "0.025"},
    {"file_name", "tmp/gr.txt"},
    {"multistate", "true"},
    {"multistate_aggregate", "false"},
  }));
  EXPECT_EQ(mc.configuration().model_params().select("cutoff").value(0), 3.);
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(1e6);
  MonteCarlo mc3 = test_serialize(mc2);
  auto pd = PairDistribution(mc3.modify(SeekModify().index("PairDistribution", mc3)[0]).modify(0));
  const grtype& gr = pd.radial(mc3.configuration());
  const Accumulator en = Energy(mc3.analyze(SeekAnalyze().index("Energy", mc3)[0])).energy();

  // compare with https://mmlapps.nist.gov/srs/LJ_PURE/md.htm
  EXPECT_EQ(gr[43].first, 1.0875);
  EXPECT_NEAR(gr[43].second[0][0], 2.65, 0.2);
  EXPECT_NEAR(en.average()/num, -5.517, 0.075);
}

}  // namespace feasst

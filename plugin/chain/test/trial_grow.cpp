#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/ideal_gas.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "steppers/include/movie.h"
#include "steppers/include/energy.h"
#include "system/include/dont_visit_model.h"
#include "chain/include/trial_grow.h"
#include "charge/include/utils.h"
#include "charge/test/charge_utils.h"

namespace feasst {

//TEST(TrialGrow, serialize) {
//  auto grow = MakeTrialGrow({
//    {{"transfer", "true"},
//     {"particle_type", "0"},
//     {"site", "0"}}});
//  auto grow2 = test_serialize(*grow);
//}

TEST(TrialGrow, angle_distribution_LONG) {
  MonteCarlo mc;
  mc.set(spce({{"physical_constants", "CODATA2010"}, {"cubic_box_length", "20"},
        {"alpha", str(5.6/20)}, {"kmax_squared", "38"}, {"dual_cut", "3.2"}}));
  mc.get_system()->get_configuration()->add_particle_of_type(0);
  mc.set(MakeThermoParams({{"beta", "1"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialGrow({
    {{"default_reference_index", "0"}, {"default_num_steps", "4"},
     {"bond", "true"}, {"particle_type", "0"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
    {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}}}));
  mc.add(MakeMovie({{"file_name", "tmp/ang"}}));
  mc.add(MakeEnergy());
  mc.attempt(1e5);
  EXPECT_NEAR(mc.analyze(mc.num_analyzers()-1).accumulator().average(), -0.087910, 0.00001);
  //INFO(mc.analyze(mc.num_analyzers()-1).accumulator().average());
}

TEST(TrialGrow, bond_harmonic) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.add(MakeConfiguration({{"cubic_box_length", "10"},
    {"particle_type0", "../plugin/chain/test/data/dimer_harmonic.fstprt"},
    {"add_particles_of_type0", "1"}}));
  mc.add(MakePotential(MakeIdealGas()));
  mc.set(MakeThermoParams({{"beta", "1"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialGrow({
    {{"bond", "true"}, {"particle_type", "0"}, {"mobile_site", "1"}, {"anchor_site", "0"}}}));
  mc.add(MakeMovie({{"file_name", "tmp/bond"}}));
  while (mc.trial(0).num_success() < 2) {
    mc.attempt(1);
  }
}

TEST(TrialGrow, file) {
  auto trial = MakeTrialGrowFile({
    {"file_name", "../plugin/chain/test/data/dimer_grow_file.txt"}});
  EXPECT_EQ(static_cast<int>(trial->trials().size()), 5);
  EXPECT_EQ(trial->trials()[0]->weight(), 50);
  EXPECT_EQ(trial->trials()[1]->weight(), 50);
  EXPECT_EQ(trial->trials()[2]->weight(), 100);
  EXPECT_EQ(trial->trials()[3]->weight(), 0.1);
  EXPECT_EQ(trial->trials()[4]->weight(), 100);
  EXPECT_EQ(trial->trials()[0]->stage(0).num_steps(), 3);
  EXPECT_EQ(trial->trials()[0]->stage(1).num_steps(), 2);
  EXPECT_EQ(trial->trials()[1]->stage(0).num_steps(), 3);
  EXPECT_EQ(trial->trials()[1]->stage(1).num_steps(), 2);
  EXPECT_EQ(trial->trials()[0]->stage(0).reference(), 2);
  EXPECT_EQ(trial->trials()[4]->stage(0).reference(), 15);
}

}  // namespace feasst

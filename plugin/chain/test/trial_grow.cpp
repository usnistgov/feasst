#include "utils/test/utils.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "steppers/include/movie.h"
#include "steppers/include/energy.h"
#include "system/include/dont_visit_model.h"
#include "chain/include/trial_grow.h"
#include "ewald/include/utils.h"

namespace feasst {

TEST(TrialGrow, serialize) {
  auto grow = MakeTrialGrow({
    {{"transfer", "true"},
     {"particle_type", "0"},
     {"site", "0"}}});
  Trial grow2 = test_serialize(*grow);
}

TEST(TrialGrow, angle_distribution_LONG) {
  MonteCarlo mc;
  mc.set(spce({{"physical_constants", "CODATA2010"}, {"cubic_box_length", "20"},
        {"alphaL", "5.6"}, {"kmax_squared", "38"}, {"dual_cut", "3.2"}}));
  mc.get_system()->get_configuration()->add_particle_of_type(0);
  mc.set(MakeThermoParams({{"beta", "1"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialGrow({
    {{"bond", "true"}, {"particle_type", "0"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
    {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}}},
    {{"reference_index", "0"}, {"num_steps", "4"}}));
  mc.add(MakeMovie({{"file_name", "tmp/ang"}}));
  mc.add(MakeEnergy());
  mc.attempt(1e5);
  EXPECT_NEAR(mc.analyze(mc.num_analyzers()-1).accumulator().average(), -0.087910, 0.00001);
  //INFO(mc.analyze(mc.num_analyzers()-1).accumulator().average());
}

}  // namespace feasst

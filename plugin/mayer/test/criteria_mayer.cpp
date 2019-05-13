#include "utils/test/utils.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_factory.h"
#include "mayer/include/criteria_mayer.h"
#include "system/test/system_test.h"

namespace feasst {

TEST(CriteriaMayer, ljb2) {
  seed_random_by_date();
  // seed_random();
  System system = default_system();
  system.add_to_reference(hs_potential());
  TrialFactory trials;
  auto translate = MakeTrialTranslate({{"max_move", "0.5"}});
  /// HWH notes: does this need a max?
  translate->set_weight(0.75);
  trials.add(translate);
  const int nTrialsEq = 1e4, nTrials = 1e4;
  //const int nTrialsEq = 1e6, nTrials = 1e6;
  Configuration * config = system.get_configuration();
  config->set_model_param("cutoff", 0, NEAR_INFINITY);
  EXPECT_EQ(config->model_params().cutoff().value(0), NEAR_INFINITY);
  const double boxl = 2*(config->model_params().cutoff().value(0));
  config->set_domain(Domain().set_cubic(boxl));
  std::cout << "boxl " << boxl << std::endl;
  CriteriaMayer criteria;
  criteria.set_beta(1.);
  criteria.add_chemical_potential(-2.775);
  criteria.set_current_energy(system.energy());
  Random random;
  for (int iTrial = 0; iTrial < nTrialsEq + nTrials; ++iTrial) {
    trials.attempt(&criteria, &system);
  }
  std::cout << "a " << criteria.second_virial() << std::endl;
  EXPECT_NEAR(-5.3, criteria.second_virial(), 15);
  EXPECT_GT(std::abs(2.0944-criteria.second_virial()), 0.0001); // HS value

  std::shared_ptr<Criteria> crit2 = test_serialize<CriteriaMayer, Criteria>(criteria);
  EXPECT_EQ(crit2->beta(), 1.);
}

}  // namespace feasst

#include <gtest/gtest.h>
#include "core/include/trial_translate.h"
#include "core/include/trial_factory.h"
#include "core/include/criteria_mayer.h"
#include "core/test/system_test.h"

namespace feasst {

TEST(CriteriaMayer, ljb2) {
  seed_random_by_date();
  // seed_random();
  System system = default_system();
  system.add_to_reference(hs_potential());
  TrialFactory trials;
  auto translate = std::make_shared<TrialTranslate>();
  translate->set_weight(0.75);
  translate->set_tunable_param_max(3);
  translate->set_max_move(0.5);
  trials.add(translate);
  const int nTrialsEq = 1e4, nTrials = 1e4;
  //const int nTrialsEq = 1e6, nTrials = 1e6;
  Configuration * config = system.get_configuration();
  config->set_model_param("cutoff", 0, NEAR_INFINITY);
  EXPECT_EQ(config->unique_types().model_params().cutoff().value(0), NEAR_INFINITY);
  const double boxl = 2*(config->unique_types().model_params().cutoff().value(0));
  config->set_domain(Domain().set_cubic(boxl));
  std::cout << "boxl " << boxl << std::endl;
  CriteriaMayer criteria;
  criteria.set_beta(1.);
  criteria.add_activity(exp(-2.775));
  criteria.set_running_energy(system.energy());
  Random random;
  for (int iTrial = 0; iTrial < nTrialsEq + nTrials; ++iTrial) {
    trials.attempt(&criteria, &system);
  }
  std::cout << "a " << criteria.second_virial() << std::endl;
  EXPECT_NEAR(-5.3, criteria.second_virial(), 15);
  EXPECT_GT(std::abs(2.0944-criteria.second_virial()), 0.0001); // HS value
}

}  // namespace feasst

#include <gtest/gtest.h>
#include "core/include/trial_translate.h"
#include "core/include/trial_factory.h"
#include "core/include/criteria_mayer.h"

TEST(CriteriaMayer, ljb2) {
  feasst::seed_random_by_date();
//  feasst::seed_random();
  feasst::System system;
  system.default_system();
//  system.model()->cut_distance = 99999999;
  feasst::TrialFactory trials;
  auto translate = std::make_shared<feasst::TrialTranslate>();
  translate->set_weight(0.75);
//  translate->set_max_move(3);
  trials.add(translate);
  const int nTrialsEq = 1e4, nTrials = 1e4;
  //const int nTrialsEq = 1e6, nTrials = 1e6;
  feasst::Configuration * config = system.configurationByPart(0);
  config->set_model_param("cutoff", 0, feasst::NEAR_INFINITY);
  EXPECT_EQ(config->unique_types().model_params().cutoff().value(0), feasst::NEAR_INFINITY);
  const double boxl = 2*(config->unique_types().model_params().cutoff().value(0));
  //const double boxl = 2*system.model()->cut_distance;
  config->set_domain(feasst::DomainCuboid().set_cubic(boxl));
  std::cout << "boxl " << boxl << std::endl;
  feasst::CriteriaMayer criteria;
  criteria.set_beta(1.);
  criteria.add_activity(exp(-2.775));
  criteria.set_running_energy(system.energy());
  feasst::Random random;

  //system.set_running_energy(system.energy());
  for (int iTrial = 0; iTrial < nTrialsEq + nTrials; ++iTrial) {
    trials.attempt(&criteria, &system);
  }
  std::cout << "a " << criteria.second_virial() << std::endl;
  EXPECT_NEAR(-5.3, criteria.second_virial(), 10);
}


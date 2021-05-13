#include <cmath>
#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/hard_sphere.h"
#include "system/include/utils.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/monte_carlo.h"
#include "mayer/include/criteria_mayer.h"

namespace feasst {

TEST(CriteriaMayer, ljb2) {
  System system = two_particle_system();
  system.add_to_reference(MakePotential(MakeHardSphere()));
  auto translate = MakeTrialTranslate({{"new_only", "true"}, {"reference_index", "0"}, {"weight", "0.75"}});
  //auto translate = MakeTrialTranslate({{"tunable_param", "0.5"}});
  /// HWH notes: does this need a max?
  const int nTrialsEq = 1e4, nTrials = 1e4;
  //const int nTrialsEq = 1e6, nTrials = 1e6;
  Configuration * config = system.get_configuration();
  config->set_model_param("cutoff", 0, NEAR_INFINITY);
  EXPECT_EQ(config->model_params().cutoff().value(0), NEAR_INFINITY);
  const double boxl = 2*(config->model_params().cutoff().value(0));
  config->set_side_lengths(Position().set_vector({boxl, boxl, boxl}));
  std::cout << "boxl " << boxl << std::endl;
  system.set(MakeThermoParams({{"beta", "1."},
    {"chemical_potential", "-2.775"}}));
  CriteriaMayer criteria;
  criteria.set_current_energy(system.energy());
  RandomMT19937 random;
  for (int iTrial = 0; iTrial < nTrialsEq + nTrials; ++iTrial) {
    translate->attempt(&criteria, &system, &random);
  }
  std::cout << "a " << criteria.second_virial() << std::endl;
  EXPECT_NEAR(-5.3, criteria.second_virial(), 15);
  EXPECT_GT(std::abs(2.0944-criteria.second_virial()), 0.0001); // HS value

  std::shared_ptr<Criteria> crit2 = test_serialize<CriteriaMayer, Criteria>(criteria);
  EXPECT_EQ(system.thermo_params().beta(), 1.);
}

CriteriaMayer ljb2(const int trials) {
  MonteCarlo mc;
  { // initialize system
    System lj = lennard_jones({{"cubic_box_length", "1000"}});
    Configuration * config = lj.get_configuration();
    config->set_model_param("cutoff", 0, config->domain().side_length(0)/2.);
    config->add_particle_of_type(0);
    config->add_particle_of_type(0);
    config->update_positions({{0, 0, 0}, {1, 1, 1}});
    mc.set(lj);
  }
  mc.add_to_reference(MakePotential(MakeHardSphere()));
  mc.set(MakeThermoParams({{"beta", "1."}, {"chemical_potential", "-2.775"}}));
  mc.set(MakeCriteriaMayer());
  mc.add(MakeTrialTranslate({{"new_only", "true"}, {"reference_index", "0"}, {"weight", "0.75"}}));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(trials);
  std::stringstream ss;
  mc2.criteria().serialize(ss);
  CriteriaMayer mayer(ss);
  return mayer;
}

TEST(MonteCarlo, ljb2) {
  CriteriaMayer mayer = ljb2(1e4);
  INFO("a " << mayer.second_virial());
  EXPECT_NEAR(-5.3, mayer.second_virial(), 20);
  EXPECT_GT(std::abs(2.0944-mayer.second_virial()), 0.0001); // HS value
}

TEST(MonteCarlo, ljb2_LONG) {
  CriteriaMayer mayer = ljb2(1e7);
  INFO("a " << mayer.second_virial());
  EXPECT_NEAR(-5.3, mayer.second_virial(), 0.2);
  EXPECT_GT(std::abs(2.0944-mayer.second_virial()), 0.0001); // HS value
}

}  // namespace feasst

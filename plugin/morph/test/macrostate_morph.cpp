#include "utils/test/utils.h"
//#include "math/include/random_mt19937.h"
#include "system/include/utils.h"
#include "system/include/long_range_corrections.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "morph/include/macrostate_morph.h"

namespace feasst {

TEST(MacrostateMorph, lj) {
  System system = lennard_jones({{"particle", "forcefield/data.lj"}});
  Configuration * config = system.get_configuration();
  config->add_particle_type(install_dir() + "/forcefield/data.lj", "0.25");
  config->set_model_param("sigma", 1, 0.25);
  config->set_model_param("cutoff", 1, 1.0);
  config->add_particle_type(install_dir() + "/forcefield/data.lj", "0.5");
  config->set_model_param("sigma", 1, 0.5);
  config->set_model_param("cutoff", 1, 1.0);
  config->add_particle_type(install_dir() + "/forcefield/data.lj", "0.75");
  config->set_model_param("sigma", 1, 0.75);
  config->set_model_param("cutoff", 1, 1.0);
  const std::vector<std::vector<int> > grow_sequence = {{1}, {2}, {3}, {0}};
  auto macro = MakeMacrostateMorph(grow_sequence,
    Histogram({{"width", str(1./grow_sequence.size())}, {"max", "10"}}));
  auto criteria = MakeMetropolis();

  Acceptance accept;
  EXPECT_EQ(macro->value(system, *criteria, accept), 0.);
  config->add_particle_of_type(1);
  EXPECT_EQ(macro->value(system, *criteria, accept), 0.25);
  TrialSelectParticle select;
  select.select_particle(0, *config);
  config->remove_particle(select.mobile());
  config->add_particle_of_type(2);
  EXPECT_EQ(macro->value(system, *criteria, accept), 0.5);
  select.select_particle(0, *config);
  config->remove_particle(select.mobile());
  config->add_particle_of_type(3);
  EXPECT_EQ(macro->value(system, *criteria, accept), 0.75);
  select.select_particle(0, *config);
  config->remove_particle(select.mobile());
  config->add_particle_of_type(0);
  EXPECT_EQ(macro->value(system, *criteria, accept), 1.);
  config->add_particle_of_type(1);
  config->update_positions({{0, 0, 0}, {1, 1, 1}});
  EXPECT_EQ(macro->value(system, *criteria, accept), 1.25);
  accept.add_to_macrostate_shift(-1);
  accept.set_macrostate_shift_type(1);
  EXPECT_EQ(macro->value(system, *criteria, accept), 1.);

  DEBUG("energy");
  DEBUG(system.energy());
  DEBUG(feasst_str(system.stored_energy_profile()));
}

}  // namespace feasst

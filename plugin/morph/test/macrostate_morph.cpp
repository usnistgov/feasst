#include "utils/test/utils.h"
//#include "math/include/random_mt19937.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/potential.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "morph/include/macrostate_morph.h"

namespace feasst {

TEST(MacrostateMorph, lj) {
  auto conf = MakeConfiguration({{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.fstprt"}});
  conf->add_particle_type(install_dir() + "/particle/lj.fstprt", "0.25");
  conf->set_model_param("sigma", 1, 0.25);
  conf->set_model_param("cutoff", 1, 1.0);
  conf->add_particle_type(install_dir() + "/particle/lj.fstprt", "0.5");
  conf->set_model_param("sigma", 1, 0.5);
  conf->set_model_param("cutoff", 1, 1.0);
  conf->add_particle_type(install_dir() + "/particle/lj.fstprt", "0.75");
  conf->set_model_param("sigma", 1, 0.75);
  conf->set_model_param("cutoff", 1, 1.0);
  System system;
  system.add(conf);
  Configuration * config = system.get_configuration();
  system.add(MakePotential(MakeLennardJones()));
  system.add(MakePotential(MakeLongRangeCorrections()));
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

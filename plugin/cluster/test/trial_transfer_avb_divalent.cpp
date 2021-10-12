#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/metropolis.h"
#include "cluster/include/trial_transfer_avb_divalent.h"
#include "cluster/include/energy_map_all.h"

namespace feasst {

// Seems there is an issue with AVBDivalent prefactors
TEST(TrialTransferAVBDivalent, add_remove) {
  System system;
  system.add(*MakeConfiguration({{"cubic_box_length", "8"},
    {"particle_type0", "../forcefield/lj.fstprt"},
    {"particle_type1", "../forcefield/atom.fstprt"}}));
  const Configuration& config = system.configuration();
  system.add(MakePotential(MakeLennardJones(),
    MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  system.energy();
  system.finalize();
  auto ran = MakeRandomMT19937();
  //ran = MakeRandomMT19937({{"seed", "1591972002"}});
  //ran = MakeRandomMT19937({{"seed", "1580154124"}});
  //ran = MakeRandomMT19937({{"seed", "1607546467"}});
  auto metropolis = MakeMetropolis();
  system.set(MakeThermoParams({
    {"beta", str(1e-6)},
    {"chemical_potential0", "2"},
    {"chemical_potential1", "2"}}));
  system.add(MakeNeighborCriteria({
    {"maximum_distance", "1.5"},
    {"minimum_distance", "1"},
    {"site_type0", "0"},
    {"site_type1", "1"},
    {"potential_index", "0"}}));
  const double vol_av = system.neighbor_criteria(0).volume(config.dimension());

  auto add = MakeTrialAddAVBDivalent({
    {"neighbor_index", "0"},
    {"particle_type", "0"},
    {"particle_type_a", "1"},
    {"particle_type_b", "1"}});
  add->precompute(metropolis.get(), &system);
  add->attempt(metropolis.get(), &system, ran.get());
  config.check();
  if (config.num_particles() != 3) return;
  EXPECT_EQ(config.particle(0).type(), 0);
  EXPECT_EQ(config.particle(1).type(), 1);
  EXPECT_EQ(config.particle(2).type(), 1);
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  DEBUG(config.particle(2).site(0).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  DEBUG(config.particle(2).site(0).position().str());
  DEBUG(config.particle(3).site(0).position().str());
  double delta = add->accept().energy_new() - add->accept().energy_old();
  EXPECT_NEAR(add->accept().ln_metropolis_prob(),
    std::log(config.domain().volume()/1.)
    +std::log(vol_av/2.)
    +std::log(vol_av/1.)
    -system.thermo_params().beta()*delta
    +system.thermo_params().beta_mu(0)
    +2*system.thermo_params().beta_mu(1),
    1e-14);

  system.energy();
  system.finalize();

  DEBUG("**begin remove test**");

  auto remove = MakeTrialRemoveAVBDivalent({
    {"neighbor_index", "0"},
    {"particle_type", "0"},
    {"particle_type_a", "1"},
    {"particle_type_b", "1"}});
  remove->precompute(metropolis.get(), &system);
  EXPECT_EQ(remove->stage(0).select().particle_type(), 0);
  EXPECT_EQ(remove->stage(1).select().particle_type(), 1);
  EXPECT_EQ(remove->stage(2).select().particle_type(), 1);
  remove->attempt(metropolis.get(), &system, ran.get());
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  DEBUG(config.particle(2).site(0).position().str());
  delta = - remove->accept().energy_old();
  EXPECT_NEAR(remove->accept().ln_metropolis_prob(),
    -std::log(config.domain().volume()/1.)
    -std::log(vol_av/2.)
    -std::log(vol_av/1.)
    -system.thermo_params().beta()*delta
    -system.thermo_params().beta_mu(0)
    -2*system.thermo_params().beta_mu(1),
    1e-14);
}

}  // namespace feasst

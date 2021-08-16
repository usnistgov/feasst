#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/dont_visit_model.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/tune.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/flat_histogram.h"
#include "ewald/include/utils.h"
#include "ewald/include/charge_screened.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/energy_map_all_criteria.h"
#include "cluster/include/energy_map_neighbor_criteria.h"
#include "chain/include/trial_grow.h"
#include "chain/include/check_rigid_bonds.h"

namespace feasst {

double energy_av44(const int macro, const MonteCarlo& mc) {
  return mc.analyzers().back()->analyzers()[macro]->accumulator().average();
}

// HWH add num steps to spce fh test for DCCB diagnosis
MonteCarlo test_spce_avb_grow_fh(std::shared_ptr<Bias> bias,
    const std::string avb_type,
    const int num_steps = 1,
    bool test = true,
    const int min = 1,
    const int max = 5,
    const int steps_per = 1e3) {
  bool avb = true;
  if (avb_type == "none") avb = false;
  DEBUG(bias->class_name());
  MonteCarlo mc;
  // mc.set(MakeRandomMT19937({{"seed", "123"}}));
  argtype spce_args = {{"physical_constants", "CODATA2010"},
                       {"cubic_box_length", "20"},
                       {"alphaL", "5.6"},
                       {"kmax_squared", "38"},
                       //{"table_size", "0"},
                      };
  int ref = -1;
  if (num_steps > 1) {
    //spce_args.insert({"dual_cut", str(10)});
    spce_args.insert({"dual_cut", str(3.16555789)});
    ref = 0;
  }
  mc.set(spce(spce_args));
  mc.get_system()->get_configuration()->add_particle_of_type(0);
  if (avb) {
    auto ncrit = MakeNeighborCriteria({{"maximum_distance", "10"}, {"minimum_distance", "3.2"}, {"site_type0", "0"}, {"site_type1", "0"}, {"potential_index", "1"}});
    mc.add(ncrit);
    auto pot = MakePotential(
      MakeLennardJones(),
      //MakeModelTwoBodyFactory({MakeLennardJones(), MakeChargeScreened({{"table_size", "0"}})}),
      //MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))//,
      MakeVisitModel(MakeVisitModelInner(MakeEnergyMapNeighborCriteria(ncrit)))//,
      //MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAllCriteria(ncrit)))//,
      //{{"table_size", "1e6"}}
    );
    mc.set(1, pot);
    mc.add(MakePotential(MakeChargeScreened({{"table_size", "0"}})));
    mc.add_to_reference(MakePotential(MakeDontVisitModel()));
    //mc.add_to_reference(pot);
  }
  const double beta = 1/kelvin2kJpermol(525, mc.configuration()); // mol/kJ
  mc.set(MakeThermoParams({{"beta", str(beta)},
     {"chemical_potential", str(-8.14/beta)}}));
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", str(max)}, {"min", str(min)}})),
    bias);
  mc.set(criteria);
//  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.275"}}));
//  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "50."}}));
  if (avb) {
    mc.add(MakeTrialGrow(
      {
        //{{"transfer", "true"}, {"particle_type", "0"}, {"weight", "4"}, {"site", "0"}},
        {{avb_type, "true"}, {"particle_type", "0"}, {"weight", "4"}, {"site", "0"}, {"neighbor_index", "0"}, {"target_particle_type", "0"}, {"target_site", "0"}},
        {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
        {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}},
      },
      {{"num_steps", "1"}, {"reference_index", "0"}}
    ));
  }
  if (avb_type != "transfer_avb") {
    mc.add(MakeTrialTransfer({
      {"particle_type", "0"},
      {"weight", "4"},
      {"reference_index", str(ref)},
      {"num_steps", str(num_steps)}}));
  }
  SeekNumParticles(min).with_thermo_params({{"beta", "1"}, {"chemical_potential", "1"}}).with_metropolis().run(&mc);
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/spce_fh"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}, {"tolerance", str(1e-6)}}));
  mc.add(MakeCheckRigidBonds({{"steps_per", str(steps_per)}}));
  mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
  mc.add(MakeCriteriaWriter({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/spce_crit.txt"}}));
  auto energy = MakeEnergy({
    {"file_name", "tmp/spce_fh_energy"},
    {"steps_per_update", "1"},
    {"steps_per_write", str(steps_per)},
    {"multistate", "true"}});
  mc.add(energy);
  //MonteCarlo mc2 = test_serialize(mc);
  mc.run_until_complete();

  if (!test) return mc;

  EXPECT_LE(mc.system().configuration().num_particles(), 5);

  // known values of lnpi and energy
  const std::vector<std::vector<double> > lnpi_srsw = {
    {-1.78421656875553, 0.015},
    {-1.47899656875553, 0.015},
    {-1.44977656875553, 0.015},
    {-1.57981656875553, 0.015},
    {-1.81051656875553, 0.015}};
//    {-2.7207, 0.015},
//    {-1.8523, 0.015},
//    {-1.54708, 0.016},
//    {-1.51786, 0.015},
//    {-1.6479, 0.015},
//    {-1.8786, 0.03}};
  const std::vector<std::vector<double> >  en_srsw = {
//    {0, 1e-13},
    {-0.0879115, 1.1293158298007674394e-05},
    {-2.326, 0.12},
    {-6.806, 0.24},
    {-13.499, 0.5},
    {-22.27, 1.0}};

  FlatHistogram fh(mc.criteria());
  const LnProbability& lnpi = fh.bias().ln_prob();
  for (int macro = 0; macro < lnpi.size(); ++macro) {
    EXPECT_NEAR(lnpi.value(macro), lnpi_srsw[macro][0],
      15*lnpi_srsw[macro][1]);
//      if (bias->class_name() == "TransitionMatrix") {
      const double en_std = std::sqrt(std::pow(en_srsw[macro][1], 2) +
        std::pow(energy->energy().block_stdev(), 2));
      EXPECT_NEAR(energy_av44(macro, mc), en_srsw[macro][0], 15.*en_std);
//      }
  }

  return mc;
}

TEST(MonteCarlo, spce_fh2_LONG) {
  //for (std::string avb_type : {"regrow_avb4"}) {
  //for (std::string avb_type : {"regrow_avb2"}) {
  //for (std::string avb_type : {"transfer_avb"}) {
  //for (std::string avb_type : {"none"}) {
  //for (std::string avb_type : {"transfer_avb", "regrow_avb2", "regrow_avb4"}) {
  for (std::string avb_type : {"none", "transfer_avb", "regrow_avb2", "regrow_avb4"}) {
    DEBUG(avb_type);
    test_spce_avb_grow_fh(MakeTransitionMatrix({{"min_sweeps", "100"}}), avb_type);
  }
}

// Seems there is an issue with TrialGrow transfer_avb spce
TEST(TrialGrow, transfer_avb_spce) {
  System system = spce({{"physical_constants", "CODATA2010"},
                        {"cubic_box_length", "20"},
                        {"alphaL", "5.6"},
                        {"kmax_squared", "38"},
                        {"add_particles_of_type0", "1"}});
  //system.get_configuration()->add_particle_of_type(0);
  auto ncrit = MakeNeighborCriteria({{"maximum_distance", "10"}, {"minimum_distance", "3.2"}, {"site_type0", "0"}, {"site_type1", "0"}, {"potential_index", "1"}});
  system.add(ncrit);
  auto pot = MakePotential(MakeLennardJones(),
    MakeVisitModel(MakeVisitModelInner(MakeEnergyMapNeighborCriteria(ncrit))));
  system.set_unoptimized(1, pot);
  system.add(MakePotential(MakeChargeScreened({{"table_size", "0"}})));
  system.add_to_reference(MakePotential(MakeDontVisitModel()));
  const Configuration& config = system.configuration();
  system.energy();
  system.finalize();
  auto ran = MakeRandomMT19937();
  //ran = MakeRandomMT19937({{"seed", "123"}});
  //ran = MakeRandomMT19937({{"seed", "1591972002"}});
  //ran = MakeRandomMT19937({{"seed", "1628624844"}});
  ran = MakeRandomMT19937({{"seed", "1628878165"}});
  auto metropolis = MakeMetropolis();
  system.set(MakeThermoParams({
    {"beta", "0.1"},
    {"chemical_potential", "1"}}));
  system.add(ncrit);
  const double vol_av = system.neighbor_criteria(0).volume(config.dimension());

  DEBUG("vol_av: " << vol_av);
  auto grow = MakeTrialGrow(
    {
      {{"transfer_avb", "true"}, {"particle_type", "0"}, {"weight", "4"}, {"site", "0"}, {"neighbor_index", "0"}, {"target_particle_type", "0"}, {"target_site", "0"}},
      {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
      {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}},
    },
    {{"num_steps", "1"}, {"reference_index", "0"}}
  );
  grow->precompute(metropolis.get(), &system);
  double en_old = metropolis->current_energy();
  bool accepted = grow->attempt(metropolis.get(), &system, 0, ran.get());
  EXPECT_EQ(grow->num(), 2);
  system.check();
  EXPECT_GT(config.num_particles(), 0);
  DEBUG(config.num_particles());
  if (config.num_particles() != 2) return;
  ASSERT(accepted, "er");
  EXPECT_EQ(config.particle(0).type(), 0);
  EXPECT_EQ(config.particle(1).type(), 0);
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(0).site(1).position().str());
  DEBUG(config.particle(0).site(2).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  DEBUG(config.particle(1).site(1).position().str());
  DEBUG(config.particle(1).site(2).position().str());
  double delta = grow->trial(0).accept().energy_new() - en_old;
  DEBUG("deltaU " << delta);
  DEBUG(grow->trial(0).accept().ln_metropolis_prob());
  EXPECT_NEAR(grow->trial(0).accept().ln_metropolis_prob(),
    std::log(vol_av/2.)
    -system.thermo_params().beta()*delta
    +system.thermo_params().beta_mu(0),
    1e-14);

//  system.energy();
//  system.finalize();
  DEBUG("en: " << system.stored_energy());
  DEBUG("en: " << metropolis->current_energy());
  en_old = metropolis->current_energy();

  DEBUG("**begin remove test**");

  accepted = grow->attempt(metropolis.get(), &system, 1, ran.get());
  //ASSERT(accepted, "er");
  delta = - grow->trial(1).accept().energy_old();
  DEBUG("deltaU: " << delta);
  EXPECT_NEAR(grow->trial(1).accept().ln_metropolis_prob(),
    -std::log(vol_av/2.)
    -system.thermo_params().beta()*delta
    -system.thermo_params().beta_mu(0),
    1e-14);

  en_old = metropolis->current_energy();
  if (accepted) return;

  DEBUG("**add a third**");
  accepted = grow->attempt(metropolis.get(), &system, 0, ran.get());
  DEBUG("old en " << en_old);
  DEBUG("new en " << grow->trial(0).accept().energy_new());
  delta = grow->trial(0).accept().energy_new() - en_old;
  DEBUG("deltaU: " << delta);
  EXPECT_DOUBLE_EQ(grow->trial(0).accept().ln_metropolis_prob(),
    std::log((2./3.)*vol_av/2.)
    -system.thermo_params().beta()*delta
    +system.thermo_params().beta_mu(0));

  INFO("** attempt avb4 **");
  if (config.num_particles() != 3) return;
  auto avb4 = MakeTrialGrow(
    {
      {{"regrow_avb4", "true"}, {"particle_type", "0"}, {"weight", "4"}, {"site", "0"}, {"neighbor_index", "0"}, {"target_particle_type", "0"}, {"target_site", "0"}},
      {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
      {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}},
    },
    {{"num_steps", "1"}, {"reference_index", "0"}}
  );
  avb4->precompute(metropolis.get(), &system);
  en_old = metropolis->current_energy();
  accepted = avb4->attempt(metropolis.get(), &system, 0, ran.get());
  INFO("accepted? " << accepted);
  INFO("new energy " << avb4->trial(0).accept().energy_new());
  delta = avb4->trial(0).accept().energy_new() - en_old;
//  EXPECT_NEAR(avb4->trial(0).accept().ln_metropolis_prob(),
//    std::log(vol_av/2.)
//    -system.thermo_params().beta()*delta
//    +system.thermo_params().beta_mu(0),
//    1e-14);
  system.check();
}

}  // namespace feasst

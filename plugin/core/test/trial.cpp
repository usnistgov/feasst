#include <gtest/gtest.h>
#include "core/include/trial_translate.h"
#include "core/include/trial_transfer.h"
#include "core/include/trial_factory.h"
#include "core/include/criteria_metropolis.h"
#include "core/include/criteria_flat_histogram.h"
#include "core/include/macrostate_num_particles.h"
#include "core/include/bias_wang_landau.h"
#include "core/include/histogram.h"
#include "core/include/utils_io.h"
#include "core/include/accumulator.h"

class TestLJSystem {
 public:
  feasst::System system;
  feasst::TrialFactory trials;
  feasst::CriteriaMetropolis criteria;
  double add_probability = 0.5;
  double box_length = 8;
  double beta = 1./1.2;
  double activ = exp(-2.775);
  feasst::Configuration * config;

  void init() {
    system.add_configuration(feasst::Configuration());
    config = system.configuration(0);
    config->add_particle_type("../forcefield/data.lj");
    auto transfer = std::make_shared<feasst::TrialTransfer>();
    transfer->set_weight(0.25);
    transfer->set_add_probability(add_probability);
    trials.add(transfer);
    auto translate = std::make_shared<feasst::TrialTranslate>();
    translate->set_weight(0.75);
    trials.add(translate);
    config->set_domain(feasst::Domain().set_cubic(box_length));
    criteria.set_beta(beta);
    criteria.add_activity(activ);
    criteria.set_running_energy(system.energy());
  }
};

TEST(Trial, MCbenchmark) {
  // feasst::seed_random_by_date();
  feasst::seed_random(1346867550);
  TestLJSystem test;
  test.add_probability = 1.;
  test.activ = 100;
  test.init();
  // const int nTrials = 1e6;  // about 3 seconds
  const int nTrials = 1e4;
  const int nMol = 50;

  // seek 50 particles
  while (test.config->num_particles() < nMol) {
    DEBUG("n " << test.config->num_particles());
    test.trials.attempt(&test.criteria, &test.system, 0);
  }

  for (int iTrial = 0; iTrial < nTrials; ++iTrial) {
    test.trials.attempt(&test.criteria, &test.system, 1);
    if (iTrial%10000==0) std::cout << "t" << iTrial << std::endl;
  }
}

TEST(Trial, NVTMC_SRSW) {
  feasst::seed_random_by_date();
  const double rho = 1e-3;
  const int nTrials = 1e3, nCheck = 1e3;
  // const int nTrials = 1e5, nCheck = 1e3;
  const int nMol = 500;
  TestLJSystem test;
  test.beta = 1./0.9;
  test.add_probability = 1.;
  test.activ = 100;
  test.box_length = pow(double(nMol)/rho, 1./3.);
  test.init();
  feasst::Accumulator pe;

  // seek particles
  while (test.config->num_particles() < nMol) {
    DEBUG("n " << test.config->num_particles());
    test.trials.attempt(&test.criteria, &test.system, 0);
  }

  for (int iTrial = 0; iTrial < 2*nTrials; ++iTrial) {
    test.trials.attempt(&test.criteria, &test.system, 1);
    if (iTrial >= nTrials) {
      pe.accumulate(test.system.visitor()->energy());
    }
    if (iTrial % nCheck == 0) {
      std::cout << "t" << iTrial
                << " pe " << pe.average()
                << " pe/n " << pe.average()/test.config->num_particles()
                << " n " << test.config->num_particles()
                << " vol " << test.config->domain().volume()
                << std::endl;
    }
  }
}

TEST(Trial, GCMCbenchmark) {
  feasst::seed_random(1346867550);
  TestLJSystem test;
  test.activ = 100;
  test.init();
  const int nTrials = 1e3;
  // test.nTrials = 1e7; // benchmark feasst0.5: ~3 seconds
  for (int iTrial = 0; iTrial < nTrials; ++iTrial) {
    test.trials.attempt(&test.criteria, &test.system);
    if (iTrial%100000==0) std::cout << "t" << iTrial << std::endl;
  }
}

TEST(Trial, MC) {
  feasst::seed_random_by_date();
  TestLJSystem test;
  const int nTrialsEq = 2e2, nTrials = 2e2;
  // const int nTrialsEq = 1e8, nTrials = 1e8;
  // const int nCheck = 1e2;
  const int nCheck = 1e4;
  const int nMol = 500;
  const double rho = 1e-3;
  test.box_length = pow(double(nMol)/rho, 1./3.);
  test.beta = 1./1.5;
  test.activ = exp(-2.775);
  test.init();
  EXPECT_EQ(1, test.config->num_particle_types());

  //system.set_running_energy(system.energy());
  int nTrans = 0, nIns = 0, nDel = 0;
  double peAcc = 0.;
  int inTrial = 0;
  for (int iTrial = 0; iTrial < nTrialsEq + nTrials; ++iTrial) {
    ++inTrial;
    peAcc += test.criteria.running_energy();
    if (iTrial == nTrialsEq) {
      peAcc = 0.;
      inTrial = 0;
    }
    if (iTrial % nCheck == 0) {
      std::cout << "n " << test.config->num_particles()
                << " pe " << test.criteria.running_energy()
                << " av " << peAcc/double(inTrial)
                << " peExpe " << -9.9165E-03*500 << std::endl;
      ASSERT(std::abs(test.system.energy() - test.criteria.running_energy()) < 1e-7, "mismatch");
    }
    test.trials.attempt(&test.criteria, &test.system);
  }
  std::cout << "nTrans " << nTrans << " nIns " << nIns << " nDel " << nDel << std::endl;
}

TEST(Trial, WLMC) {
  // feasst::seed_random();
  feasst::seed_random_by_date();
  feasst::System system;
  system.add_configuration(feasst::Configuration());
  feasst::Configuration * config = system.configuration(0);
  config->add_particle_type("../forcefield/data.lj");

  // trials
  feasst::TrialFactory trials;
  auto transfer = std::make_shared<feasst::TrialTransfer>();
  transfer->set_weight(0.25);
  trials.add(transfer);
  //transfer.set_add_probability(1.);
  auto translate = std::make_shared<feasst::TrialTranslate>();
  translate->set_weight(1);
  trials.add(translate);

  const int nTrialsEq = 2e2, nTrials = 2e2;
  //const int nTrialsEq = 1e8, nTrials = 1e8;
  //const int nCheck = 1e6;
  const int nCheck = 1e2;
  const double boxl = 8;
  //const double rho = 1e-3, boxl = pow(double(nMol)/rho, 1./3.);
  config->set_domain(feasst::Domain().set_cubic(boxl));
  std::cout << "boxl " << boxl << std::endl;
  feasst::CriteriaFlatHistogram criteria;
  criteria.set_beta(1./1.5);
  criteria.add_activity(exp(-2.775));
  //criteria.add_activity(exp(-2.775));
  feasst::Histogram histogram;
  const int nMol = 5;
  histogram.set_width_center(1., 0.);
  for (int i = 0; i <= nMol; ++i) {
    histogram.add(i);
  }
  auto macrostate = std::make_shared<feasst::MacrostateNumParticles>();
  macrostate->set_histogram(histogram);
  criteria.set_macrostate(macrostate);
  auto bias = std::make_shared<feasst::BiasWangLandau>();
  bias->resize(histogram);
  criteria.set_bias(bias);

  criteria.set_running_energy(system.energy());
  double peAcc = 0.;
  int inTrial = 0;
  for (int iTrial = 0; iTrial < nTrialsEq + nTrials; ++iTrial) {
    ++inTrial;
    peAcc += criteria.running_energy();
    if (iTrial == nTrialsEq) {
      peAcc = 0.;
      inTrial = 0;
    }
    if (iTrial % nCheck == 0) {
      std::cout << "n " << config->num_particles() << " pe " << criteria.running_energy() << " av " << peAcc/double(inTrial) << " peExpe " << -9.9165E-03*500 << std::endl;
      ASSERT(std::abs(system.energy() - criteria.running_energy()) < 1e-7, "mismatch");
    }
    trials.attempt(&criteria, &system);
  }
  std::cout << "nTrans " << translate->num_attempts() << " nTransfer " << transfer->num_attempts() << std::endl;
  std::cout << feasst::str(bias->ln_macro_prob()) << std::endl;
}

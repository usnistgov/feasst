#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "configuration/include/utils.h"
#include "configuration/include/domain.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/select.h"
#include "system/include/visit_model_cell.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

/// Use a 5 angstrom cut off
TEST(VisitModelCell, lj_reference_config) {
  Configuration config = lj_sample4();
  const int rcut = 2;
  for (int site_index = 0; site_index < config.num_site_types(); ++site_index) {
    config.set_model_param("cutoff", site_index, rcut);
  }
  config.add_particle_of_type(0);
  Select select(config.newest_particle_index(), config.newest_particle());
  config.remove_particle(select);
  auto domain = std::make_shared<Domain>(config.domain());
  config.set(domain);
  config.check();
  LennardJones model;
  auto cell_visit = MakeVisitModelCell({{"min_length", feasst::str(rcut)}});
  VisitModel visit;
  cell_visit->precompute(&config);
  visit.precompute(&config);
  model.compute(&config, &visit);
  model.compute(&config, cell_visit.get());
  EXPECT_NEAR(visit.energy(), cell_visit->energy(), 5e-12);
  EXPECT_NEAR(-15.076312312129405, visit.energy(), NEAR_ZERO);

  /// test energy of selection
  auto tsel = MakeTrialSelectParticle({{"particle_type", "0"}});
  RandomMT19937 random;
  tsel->random_particle(config, &select, &random);
  model.compute(select, &config, &visit);
  model.compute(select, &config, cell_visit.get());
  EXPECT_NEAR(visit.energy(), cell_visit->energy(), 5e-15);

  auto cell_visit2 = test_serialize<VisitModelCell, VisitModel>(*cell_visit);
  model.compute(select, &config, cell_visit2.get());
  EXPECT_NEAR(visit.energy(), cell_visit2->energy(), 5e-15);

  visit.check_energy(&model, &config);
  cell_visit->check_energy(&model, &config);
}

/// Use a 5 angstrom cut off
TEST(VisitModelCell, spce_reference_config) {
  Configuration config = spce_sample1();
  const int rcut = 5;
  for (int site_index = 0; site_index < config.num_site_types(); ++site_index) {
    config.set_model_param("cutoff", site_index, rcut);
  }
  config.add_particle_of_type(0);
  Select select(config.newest_particle_index(), config.newest_particle());
  config.remove_particle(select);
  auto domain = std::make_shared<Domain>(config.domain());
  config.set(domain);
  config.check();
  LennardJones model;
  auto cell_visit = MakeVisitModelCell({{"min_length", feasst::str(rcut)}});
  VisitModel visit;
  cell_visit->precompute(&config);
  visit.precompute(&config);
  model.compute(&config, &visit);
  model.compute(&config, cell_visit.get());
  EXPECT_NEAR(visit.energy(), cell_visit->energy(), 5e-12);
  EXPECT_NEAR(896.85497602741475, visit.energy(), NEAR_ZERO);

  /// test energy of selection
  auto tsel = MakeTrialSelectParticle({{"particle_type", "0"}});
  RandomMT19937 random;
  tsel->random_particle(config, &select, &random);
  model.compute(select, &config, &visit);
  model.compute(select, &config, cell_visit.get());
  EXPECT_NEAR(visit.energy(), cell_visit->energy(), 5e-14);

  visit.check_energy(&model, &config);
  cell_visit->check_energy(&model, &config);
}

//// add individual particles until reaching the point where the visitors are
//// inconsistent
//TEST(VisitModelCell, spce_reference_config_buildup) {
//  Configuration config;
//  config.add_particle_type("../forcefield/data.spce");
//  const int rcut = 5;
//  for (int site_index = 0; site_index < config.num_site_types(); ++site_index) {
//    config.set_model_param("cutoff", site_index, rcut);
//  }
//  for (int part = 0; part < 100; ++part) {
//    config.add_particle_of_type(0);
//  }
//  FileXYZ().load("../plugin/system/test/data/spce_sample_config_periodic1.xyz", &config);
//  auto domain = std::make_shared<Domain>(*config.domain());
//  domain->init_cells(rcut);
//  config.set(domain);
//  LennardJones model;
//  VisitModelCell cell_visit;
//  VisitModel visit;
//  cell_visit.precompute(&config);
//  visit.precompute(&config);
//
//  System sys;
//  sys.add(config);
//  sys.add(config);
//  Configuration * config1 = sys.get_configuration(0);
//  Configuration * config2 = sys.get_configuration(1);
//
//  // remove particles in config2
//  EXPECT_EQ(100, config1->num_particles());
//  EXPECT_EQ(100, config2->num_particles());
//  config2->remove_particles(config2->group_select(0));
//  EXPECT_EQ(0, config2->num_particles());
//
//  PerturbConfigs perturb;
//  int transfers = 0;
//  while (config1->num_particles() > 90) {
//    perturb.transfer_particle(0, &sys, 0, 1);
//    ++transfers;
//    DEBUG("transfers " << transfers);
//    DEBUG("cell list " << config2->domain()->cells(0).str());
//    EXPECT_EQ(100 - transfers, config1->num_particles());
//    EXPECT_EQ(transfers, config2->num_particles());
//    cell_visit.check_energy(model, config2);
//  }
//}

}  // namespace feasst

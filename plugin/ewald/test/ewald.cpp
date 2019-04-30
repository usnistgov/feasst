#include <gtest/gtest.h>
#include "ewald/include/ewald.h"
#include "core/test/configuration_test.h"
#include "ewald/include/ewald.h"
#include "ewald/include/model_charge_self.h"
#include "ewald/include/model_charge_intra.h"
#include "ewald/include/model_charge_screened.h"
#include "core/include/model_lj.h"
#include "core/include/long_range_corrections.h"
#include "core/include/visit_model_intra.h"
#include "core/include/system.h"
#include "core/include/perturb_translate.h"
#include "core/include/monte_carlo.h"
#include "core/include/trial_translate.h"
#include "core/include/criteria_metropolis.h"
#include "ewald/test/system_example.h"

namespace feasst {

TEST(Ewald, ewald) {
  Configuration config = spce_sample();
  config.add_model_param("alpha", 5.6/config.domain().min_side_length());
  Ewald ewald;
  ewald.set_kmax(27); // if this is changed, also update eik**0_index
  ewald.update_wave_vectors(config);
  ewald.init_wave_vector_storage(&config);
//  INFO(config.particle(0).site(0).properties().property_name()));
//  INFO(feasst::str(config.particle(0).site(0).properties().property_name()));
//  INFO(feasst::str(config.particle(0).site(0).properties().property_value()));

  ModelChargeSelf model;  // any place holder model is fine because its not used
  model.compute(&config, &ewald);
  // ewald.update_eik(config.selection_of_all(), &config);

  const std::vector<double> eik = config.particle(0).site(0).properties().property_value();
  EXPECT_NEAR(eik[0], 1, NEAR_ZERO);
  EXPECT_NEAR(eik[1], -0.069470287276879206, NEAR_ZERO);
  EXPECT_NEAR(eik[2], -0.99034775837133582, NEAR_ZERO);
  const int ry0_index = 13;
  EXPECT_NEAR(eik[ry0_index + 0], 1, NEAR_ZERO);
  EXPECT_NEAR(eik[ry0_index + 1], -0.87389397051446949, NEAR_ZERO);
  EXPECT_NEAR(eik[ry0_index - 1], eik[ry0_index + 1], NEAR_ZERO);
  const int ix0_index = 33;
  EXPECT_NEAR(eik[ix0_index + 0], 0, NEAR_ZERO);
  EXPECT_NEAR(eik[ix0_index + 1], -0.99758402111584965, NEAR_ZERO);
  EXPECT_NEAR(eik[ix0_index + 2], 0.13860489705948481, NEAR_ZERO);
  const int iz0_index = 59;
  EXPECT_NEAR(eik[iz0_index + 0], 0, NEAR_ZERO);
  EXPECT_NEAR(eik[iz0_index + 1], -0.52837486359383823, NEAR_ZERO);
  EXPECT_NEAR(eik[iz0_index - 1], 0.52837486359383823, NEAR_ZERO);

  EXPECT_NEAR(ewald.struct_fact_real()[0], -1.829963812936731, 5e-15);
  EXPECT_NEAR(ewald.struct_fact_imag()[0], 2.3263016099862206, 5e-15);

  EXPECT_NEAR(ewald.energy(), 52.1324574151071, 1e-12);

  // serialize
  std::stringstream ss;
  ewald.serialize(ss);
  auto ewald2 = Ewald().deserialize(ss);
  model.compute(&config, ewald2.get());
  EXPECT_NEAR(ewald2->energy(), 52.1324574151071, 1e-12);
}

TEST(Ewald, system) {
  System sys = spce();
  EXPECT_NEAR(-4062.47263092246, sys.energy(), 1e-10);
  EXPECT_NEAR(-3819.24971214984,
    sys.unoptimized().potentials()[0].stored_energy(), 1e-10);
  EXPECT_NEAR(23363.573774608,
    sys.unoptimized().potentials()[1].stored_energy(), 1e-10);
  EXPECT_NEAR(-6.84874714555147,
    sys.unoptimized().potentials()[2].stored_energy(), 1e-13);
  EXPECT_NEAR(-23652.08040365018,
    sys.unoptimized().potentials()[3].stored_energy(), 1e-12);
  EXPECT_NEAR(52.1324574151071,
    sys.unoptimized().potentials()[4].stored_energy(), 1e-12);
}

TEST(Ewald, revert) {
  seed_random_by_date();
  System sys = spce();
  const std::vector<double> disp = {1.43, -2.5, 0.03};
  Position trajectory(disp);
  EXPECT_NEAR(-4062.47263092246, sys.energy(), 1e-10);
  PerturbTranslate perturb;
  perturb.select_random_particle(0, sys.configuration());
  perturb.translate_selection(trajectory, &sys);
  EXPECT_GT(std::abs(-4062.47263092246 - sys.energy()), 1e-10);
  perturb.revert();
  EXPECT_NEAR(-4062.47263092246, sys.energy(), 1e-10);
}

TEST(Ewald, MC) {
  MonteCarlo mc;
  mc.set(spce());
  auto criteria = std::make_shared<CriteriaMetropolis>();
  criteria->set_beta(1.);
  mc.set(criteria);
  auto translate = MakeTrialTranslate({{"max_move", "0.1"}});
  mc.add(translate);
//  mc.attempt(100);
//  INFO(mc.get_system()->energy());
}

}  // namespace feasst

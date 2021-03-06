#include <cmath>
#include "utils/test/utils.h"
#include "configuration/include/utils.h"
#include "math/include/random_mt19937.h"
#include "system/include/model_empty.h"
#include "system/include/system.h"
#include "ewald/include/ewald.h"
#include "ewald/include/utils.h"
#include "ewald/test/system_example.h"
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

TEST(Ewald, ewald) {
  Configuration config = spce_sample1();
  auto ewald = MakeEwald({
    {"alpha", str(5.6/config.domain().min_side_length())},
    {"kmax_squared", "27"}
  });
  ewald->precompute(&config);

  ModelEmpty model;  // any place holder model is fine because its not used
  model.compute(&config, ewald.get());
  ewald->finalize(Select());
  // ewald.update_eik(config.selection_of_all(), &config);

  const std::vector<double> eik = config.particle(0).site(0).properties().values();
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

  EXPECT_NEAR(ewald->struct_fact_real()[0], -1.829963812936731, 5e-15);
  EXPECT_NEAR(ewald->struct_fact_imag()[0], 2.3263016099862206, 5e-15);

  const double en = 52.13245734204947;
  EXPECT_NEAR(ewald->energy(), en, 1e-12);

  EXPECT_NEAR(ewald->net_charge(config), 0., NEAR_ZERO);

  // serialize
  auto ewald2 = test_serialize<Ewald, VisitModel>(*ewald);
  model.compute(&config, ewald2.get());
  ewald2->finalize(Select());
  EXPECT_NEAR(ewald2->energy(), en, 1e-12);
}

TEST(Ewald, system) {
  const double en_lrc = -6.84874714555147;
  {
    System system = spce({{"physical_constants", "CODATA2010"},
      {"kmax_squared", "27"},
      {"xyz_file", "../plugin/configuration/test/data/spce_sample_config_periodic1.xyz"}});
    EXPECT_NEAR(-4062.472630678883, system.energy(), 1e-10);
    EXPECT_NEAR(52.1324574151071, system.potential(0).stored_energy(), 1e-12);
    EXPECT_NEAR(-3819.2497119062596, system.potential(1).stored_energy(), 1e-10);
    EXPECT_NEAR(23363.573774608, system.potential(2).stored_energy(), 1e-10);
    EXPECT_NEAR(-23652.08040365018, system.potential(3).stored_energy(), 1e-12);
    EXPECT_NEAR(en_lrc, system.potential(4).stored_energy(), 1e-13);
  }
  {
    System system = spce({{"physical_constants", "CODATA2018"},
      {"kmax_squared", "27"},
      {"xyz_file", "../plugin/configuration/test/data/spce_sample_config_periodic1.xyz"}});
    EXPECT_NEAR(-4062.4726238355574, system.energy(), 1e-10);
    EXPECT_NEAR(52.13245734204947, system.potential(0).stored_energy(), 1e-12);
    EXPECT_NEAR(-3819.2497053941993, system.potential(1).stored_energy(), 1e-10);
    EXPECT_NEAR(23363.573741866534, system.potential(2).stored_energy(), 1e-10);
    EXPECT_NEAR(-23652.080370504391, system.potential(3).stored_energy(), 1e-12);
    EXPECT_NEAR(en_lrc, system.potential(4).stored_energy(), 1e-13);
  }
}

TEST(Ewald, revert) {
  System system = spce({{"physical_constants", "CODATA2018"},
    {"kmax_squared", "27"},
    {"xyz_file", "../plugin/configuration/test/data/spce_sample_config_periodic1.xyz"}});
  const double en = -4062.4726238355574;
  EXPECT_NEAR(en, system.energy(), 1e-10);
  PerturbTranslate perturb;
  TrialSelectParticle tsel;
  RandomMT19937 random;
  tsel.sel(&system, &random);
  perturb.perturb(&system, &tsel, &random);
  EXPECT_GT(std::abs(en - system.energy()), 1e-10);
  perturb.revert(&system);
  EXPECT_NEAR(en, system.energy(), 1e-10);
}

TEST(Ewald, change_volume) {
  System system = spce();
  try {
    system.change_volume(1);
    CATCH_PHRASE("not implemented");
  }
}

}  // namespace feasst

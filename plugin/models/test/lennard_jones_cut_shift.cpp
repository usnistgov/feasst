#include <cmath>  // pow
#include "utils/test/utils.h"
#include "models/include/lennard_jones_cut_shift.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(LennardJonesCutShift, analytical) {
  Configuration config;
  config.add_particle_type("../forcefield/lj.fstprt");
  auto shift = std::make_shared<LennardJonesCutShift>();
  shift->precompute(config.model_params());
  EXPECT_NEAR(0., shift->energy(3*3, 0, 0, config.model_params()), NEAR_ZERO);
  EXPECT_NEAR(-0.010837449391761200, shift->energy(2.5*2.5, 0, 0, config.model_params()), NEAR_ZERO);
}

TEST(LennardJonesCutShift, analytical_delta) {
  Configuration config;
  config.add(MakeDeltaSigma());
  config.add_particle_type("../plugin/models/forcefield/ljdelta.fstprt");
  //auto shift = std::make_shared<LennardJonesAlpha>();
  auto shift = std::make_shared<LennardJonesCutShift>();
  shift->precompute(config.model_params());
  EXPECT_NEAR(0, shift->energy(2.9999999999999999999999*2.9999999999999999999999, 0, 0, config.model_params()), NEAR_ZERO);
  EXPECT_NEAR(-0.15885125916246515 - -0.0021747803916549908,
    shift->energy(1.2*1.2, 0, 0, config.model_params()), NEAR_ZERO);
}

TEST(LennardJonesCutShift, wca) {
  Configuration config;
  config.add_particle_type("../forcefield/lj.fstprt");
  auto wca = std::make_shared<LennardJonesCutShift>();
  ModelParams wca_params = deep_copy(config.model_params());
  EXPECT_NEAR(3., config.model_params().select("cutoff").mixed_values()[0][0], NEAR_ZERO);
  wca->set_wca(0, 0, &wca_params);
  EXPECT_NEAR(3., config.model_params().select("cutoff").mixed_values()[0][0], NEAR_ZERO);
  EXPECT_NEAR(1.122462048309373, wca_params.select("cutoff").mixed_values()[0][0], NEAR_ZERO);
  wca->precompute(wca_params);
  const double r_wca = std::pow(2, 1./6.);
  EXPECT_NEAR(0., wca->energy(r_wca*r_wca, 0, 0, wca_params), NEAR_ZERO);
  EXPECT_NEAR(1., wca->energy(1., 0, 0, wca_params), NEAR_ZERO);
}

TEST(LennardJonesCutShift, serialize) {
  Configuration config;
  config.add_particle_type("../forcefield/spce.fstprt");
  auto shift = MakeLennardJonesCutShift({{"alpha", "12"},
                                         {"hard_sphere_threshold", "0.3"}});
  shift->precompute(config.model_params());
  std::shared_ptr<Model> model2 = test_serialize<LennardJonesCutShift, Model>(*shift,
    "LennardJonesCutShift 2094 0 1 2 3 763 0.089999999999999997 713 12 -1 0 -1 1.0594630943592953 644 1 energy_at_cutoff 795 2 0 0 2 2 -2.6332331818264547e-06 -0 2 -0 -0 2 2 1 1 2 1 1 ");
}

TEST(LennardJonesCutShift, analytical_lambda) {
  Configuration config;
  config.add(MakeLambda());
  config.add(MakeDeltaSigma());
  config.add_particle_type("../plugin/models/forcefield/ljlambda.fstprt");
  auto shift = MakeLennardJonesCutShift({{"lambda", "true"}});
  shift->precompute(config.model_params());
  EXPECT_NEAR(-0.002739720872119390 - -0.001087390195827500, shift->energy(2.5*2.5, 0, 0, config.model_params()), NEAR_ZERO);
  EXPECT_NEAR(NEAR_ZERO, shift->energy(3*3, 0, 0, config.model_params()), NEAR_ZERO);
}

}  // namespace feasst

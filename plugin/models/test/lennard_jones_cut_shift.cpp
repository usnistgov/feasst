#include <cmath>  // pow
#include "utils/test/utils.h"
#include "models/include/lennard_jones_cut_shift.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(LennardJonesCutShift, analytical) {
  Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  auto shift = std::make_shared<LennardJonesCutShift>();
  shift->precompute(config.model_params());
  EXPECT_NEAR(0., shift->energy(3*3, 0, 0, config.model_params()), NEAR_ZERO);
  EXPECT_NEAR(-0.010837449391761200, shift->energy(2.5*2.5, 0, 0, config.model_params()), NEAR_ZERO);
}

TEST(LennardJonesCutShift, wca) {
  Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  auto wca = std::make_shared<LennardJonesCutShift>();
  ModelParams wca_params = config.model_params();
  EXPECT_NEAR(3., config.model_params().mixed_cutoff()[0][0], NEAR_ZERO);
  wca->set_wca(0, 0, &wca_params);
  EXPECT_NEAR(3., config.model_params().mixed_cutoff()[0][0], NEAR_ZERO);
  wca->precompute(wca_params);
  const double r_wca = std::pow(2, 1./6.);
  EXPECT_NEAR(0., wca->energy(r_wca*r_wca, 0, 0, wca_params), NEAR_ZERO);
  EXPECT_NEAR(1., wca->energy(1., 0, 0, wca_params), NEAR_ZERO);
}

TEST(LennardJonesCutShift, serialize) {
  Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  auto shift = MakeLennardJonesCutShift({{"alpha", "12"},
                                         {"hard_sphere_threshold", "0.3"}});
  shift->precompute(config.model_params());
  std::shared_ptr<Model> model2 = test_serialize<LennardJonesCutShift, Model>(*shift,
    "LennardJonesCutShift 763 0.089999999999999997 713 12 644 1 ModelParam 795 2 0 0 2 2 -2.6332331818264547e-06 -0 2 -0 -0 2 2 1 1 2 1 1 ");
}

}  // namespace feasst

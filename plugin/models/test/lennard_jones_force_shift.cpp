#include "utils/test/utils.h"
#include "models/include/lennard_jones_force_shift.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(LennardJonesForceShift, analytical) {
  auto config = MakeConfiguration({{"particle_type0", "../forcefield/lj.fstprt"}});
  EXPECT_NEAR(3., config->model_params().select("cutoff").mixed_values()[0][0], NEAR_ZERO);
  auto force_shift = std::make_shared<LennardJonesForceShift>();
  force_shift->precompute(config->model_params());
  EXPECT_NEAR(-0, force_shift->energy(3*3, 0, 0, config->model_params()), NEAR_ZERO);
  EXPECT_NEAR(-0.005365534353215080, force_shift->energy(2.5*2.5, 0, 0, config->model_params()), NEAR_ZERO);
}

TEST(LennardJonesForceShift, analytical_delta) {
  auto config = MakeConfiguration({{"particle_type0", "../plugin/models/forcefield/ljdelta.fstprt"}});
  EXPECT_NEAR(3., config->model_params().select("cutoff").mixed_values()[0][0], NEAR_ZERO);
  auto force_shift = std::make_shared<LennardJonesForceShift>();
  force_shift->precompute(config->model_params());
  EXPECT_NEAR(-0, force_shift->energy(3*3, 0, 0, config->model_params()), NEAR_ZERO);
  EXPECT_NEAR(-0.0054794417442387755 - -0.0021747803916549908 - -0.5*0.003726165748658701, force_shift->energy(2.5*2.5, 0, 0, config->model_params()), NEAR_ZERO);
}

TEST(LennardJonesForceShift, serialize) {
  auto config = MakeConfiguration({{"particle_type0", "../forcefield/spce.fstprt"}});
  auto shift = MakeLennardJonesForceShift({{"alpha", "12"},
      {"hard_sphere_threshold", "0.3"}});
  shift->precompute(config->model_params());
  std::shared_ptr<Model> model2 = test_serialize<LennardJonesForceShift, Model>(*shift,
    "LennardJonesForceShift 2094 1 0 2 3 763 0.089999999999999997 713 12 -1 0 -1 1.0594630943592953 923 energy_at_cutoff 4795 2 0 0 2 2 -2.6332331818264547e-06 -0 2 -0 -0 2 2 1 1 2 1 1 energy_deriv_at_cutoff 4795 2 0 0 2 2 3.1598766187506022e-06 0 2 0 0 2 2 1 1 2 1 1 1 ");
}

TEST(LennardJonesForceShift, analytical_lambda) {
  auto config = MakeConfiguration({{"particle_type0", "../plugin/models/forcefield/ljlambda.fstprt"}});
  auto model = MakeLennardJonesForceShift({{"lambda", "true"}});
  model->precompute(config->model_params());
  EXPECT_NEAR(-0.002739720872119390 - -0.001087390195827500 - -0.5*0.5*0.003726165748658701, model->energy(2.5*2.5, 0, 0, config->model_params()), NEAR_ZERO);
  EXPECT_NEAR(NEAR_ZERO, model->energy(3*3, 0, 0, config->model_params()), NEAR_ZERO);

}
}  // namespace feasst

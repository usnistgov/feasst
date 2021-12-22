#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "system/include/lennard_jones.h"
#include "models/include/mie.h"

namespace feasst {

TEST(Mie, analytical) {
  Configuration config;
  config.add_particle_type("../forcefield/lj.fstprt");
  auto model1 = std::make_shared<LennardJones>();
  auto model2 = std::make_shared<Mie>();
  model1->precompute(config.model_params());
  model2->precompute(config.model_params());
  EXPECT_NEAR(model1->energy(3.*3., 0, 0, config.model_params()),
              model2->energy(3.*3., 0, 0, config.model_params()), NEAR_ZERO);

  auto model3 = MakeMie({{"n", "14"}, {"m", "8"}});
  model3->precompute(config.model_params());
  std::shared_ptr<Model> model4 = test_serialize<Mie, Model>(*model3, "Mie 2094 1 0 2 -1 2905 14 8 4.9207071226910948 ");
  INFO(model4->energy(1.5*1.5, 0, 0, config.model_params()));
  EXPECT_NEAR(-0.17514250679168769, model4->energy(1.5*1.5, 0, 0, config.model_params()), NEAR_ZERO);
}

}  // namespace feasst

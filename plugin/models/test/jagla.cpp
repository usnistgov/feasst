#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "configuration/include/model_params.h"
#include "models/include/jagla.h"

namespace feasst {

TEST(Jagla, serialize) {
  auto example = MakeJagla();
  auto example_serialized = test_serialize(*example);
}

void test_energy_jagla(const double distance, const double energy,
    std::shared_ptr<Model> model, std::shared_ptr<Configuration> config) {
  EXPECT_NEAR(energy,
              model->energy(distance*distance, 0, 0, config->model_params()),
              NEAR_ZERO);
}

TEST(Jagla, analytical) {
  auto config = MakeConfiguration({{"particle_type0", "../plugin/example/particle/jagla.fstprt"}});
  auto jagla = std::make_shared<Jagla>();
  jagla->precompute(config->model_params());
  EXPECT_NEAR(1.5, config->model_params().select("lambda").value(0), NEAR_ZERO);
  EXPECT_NEAR(1, config->model_params().select("gamma").value(0), NEAR_ZERO);
  test_energy_jagla(0.999999, NEAR_INFINITY, jagla, config);
  test_energy_jagla(1, 1, jagla, config);
  test_energy_jagla(1.25, 0., jagla, config);
  test_energy_jagla(1.5, -1, jagla, config);
  test_energy_jagla(2.25, -0.5, jagla, config);
  test_energy_jagla(3., 0, jagla, config);
}

}  // namespace feasst

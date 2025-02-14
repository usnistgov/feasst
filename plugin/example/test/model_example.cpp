#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "configuration/include/model_params.h"
#include "example/include/model_example.h"

namespace feasst {

// Test that the model is serializable and remembers custom parameters.
TEST(ModelExample, serialize) {
  auto example = std::make_shared<ModelExample>(argtype({{"num_discretized_steps", "10"}}));
  auto example_serialized = test_serialize(*example);
  EXPECT_EQ(10, example_serialized.num_discretized_steps());
}

void test_energy(const double distance, const double energy,
    std::shared_ptr<Model> model, std::shared_ptr<Configuration> config) {
  EXPECT_NEAR(energy,
              model->energy(distance*distance, 0, 0, config->model_params()),
              NEAR_ZERO);
}

TEST(ModelExample, analytical) {
  auto config = MakeConfiguration({{"particle_type0", "../plugin/example/particle/jagla.fstprt"}});
  auto jagla = std::make_shared<ModelExample>();
  jagla->precompute(config->model_params());

  // Are lambda and gamma read from the fstprt file and stored appropriately?
  EXPECT_NEAR(1.5, config->model_params().select("lambda").value(0), NEAR_ZERO);
  EXPECT_NEAR(1, config->model_params().select("gamma").value(0), NEAR_ZERO);

  // Test a range of distances
  test_energy(0.999999, NEAR_INFINITY, jagla, config);
  test_energy(1, 1, jagla, config);
  test_energy(1.25, 0., jagla, config);
  test_energy(1.5, -1, jagla, config);
  test_energy(2.25, -0.5, jagla, config);
  test_energy(3., 0, jagla, config);
}

}  // namespace feasst

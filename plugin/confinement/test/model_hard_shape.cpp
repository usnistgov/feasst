#include "utils/test/utils.h"
#include "confinement/include/model_hard_shape.h"
#include "confinement/include/half_space.h"
#include "configuration/include/configuration.h"
#include "system/include/select_list.h"

namespace feasst {

TEST(ModelHardShape, half_space) {
  HalfSpace half_space({
    {"dimension", "2"},
    {"intersection", "1."},
    {"direction", "1"},
  });
  ModelHardShape model(std::make_shared<HalfSpace>(half_space));
  auto config = MakeConfiguration(MakeDomain({{"cubic_box_length", "8"}}),
    {{"particle_type", "../forcefield/data.atom"}});
  config->add_particle_of_type(0);
  const ModelParams model_params = config->model_params();
  EXPECT_LT(1e100, model.energy(config->particle(0).site(0), config.get(), model_params));
  Position pos;
  pos.set_vector({-24.23, 35.45, 1.5000001});
  config->displace_particles(config->selection_of_all(), pos);
  EXPECT_NEAR(0, model.energy(config->particle(0).site(0), config.get(), model_params), 1e-15);
  pos.set_vector({0, 0, -.00001});
  config->displace_particle(config->selection_of_all(), pos);
  EXPECT_LT(1e100, model.energy(config->particle(0).site(0), config.get(), model_params));

  std::shared_ptr<Model> model2 = test_serialize<ModelHardShape, Model>(model);
}

}  // namespace feasst

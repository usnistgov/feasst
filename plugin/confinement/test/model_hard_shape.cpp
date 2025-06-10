#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "shape/include/half_space.h"
#include "configuration/include/model_params.h"
#include "configuration/include/configuration.h"
#include "configuration/include/select.h"
#include "configuration/include/domain.h"
#include "confinement/include/model_hard_shape.h"

namespace feasst {

TEST(ModelHardShape, half_space) {
  HalfSpace half_space({
    {"dimension", "2"},
    {"intersection", "1."},
    {"direction", "1"},
  });
  ModelHardShape model(std::make_shared<HalfSpace>(half_space));
  auto config = MakeConfiguration({{"cubic_side_length", "8"},
    {"particle_type", "../particle/atom.txt"},
    {"add_particles_of_type0", "1"}});
  const ModelParams model_params = config->model_params();
  model.precompute(*config);
  EXPECT_LT(1e100, model.energy_no_wrap(config->particle(0).site(0), *config, model_params));
  Position pos;
  pos.set_vector({-24.23, 35.45, 1.5000001});
  config->displace_particles(config->selection_of_all(), pos);
  EXPECT_NEAR(0, model.energy_no_wrap(config->particle(0).site(0), *config, model_params), 1e-15);
  pos.set_vector({0, 0, -.00001});
  config->displace_particle(config->selection_of_all(), pos);
  EXPECT_LT(1e100, model.energy_no_wrap(config->particle(0).site(0), *config, model_params));

  std::shared_ptr<Model> model2 = test_serialize<ModelHardShape, Model>(model);
}

}  // namespace feasst

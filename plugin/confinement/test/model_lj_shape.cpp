#include "utils/test/utils.h"
#include "shape/include/half_space.h"
#include "configuration/include/configuration.h"
#include "configuration/include/select.h"
#include "configuration/include/domain.h"
#include "confinement/include/model_lj_shape.h"

namespace feasst {

TEST(ModelLJShape, half_space) {
  HalfSpace half_space({
    {"dimension", "2"},
    {"intersection", "1."},
    {"direction", "1"},
  });
  ModelLJShape model(std::make_shared<HalfSpace>(half_space));
  auto config = MakeConfiguration({{"cubic_side_length", "8"},
    {"particle_type", "../forcefield/atom.fstprt"},
    {"add_particles_of_type0", "1"}});
  const ModelParams model_params = config->model_params();
  model.precompute(model_params);
  std::shared_ptr<Model> model2 = test_serialize<ModelLJShape, Model>(model);

  const Site& site = config->particle(0).site(0);
  EXPECT_LT(1e100, model.energy(site.position(), site, *config, model_params));
  Position pos;
  pos.set_vector({-24.23, 35.45, 1.5000001});
  config->displace_particles(config->selection_of_all(), pos);
  EXPECT_NEAR(7.96295816296489, model.energy(site.position(), site, *config, model_params), 2e-14);
  pos.set_vector({0, 0, 2.5});
  config->displace_particle(config->selection_of_all(), pos);
  EXPECT_NEAR(0, model.energy(site.position(), site, *config, model_params), 1e-15);
  pos.set_vector({0, 0, -0.00001});
  config->displace_particle(config->selection_of_all(), pos);
  EXPECT_NEAR(3.66669086690474E-07, model.energy(site.position(), site, *config, model_params), 1e-15);

  // Without shift
  ModelLJShape model3(std::make_shared<HalfSpace>(half_space), {{"disable_shift", "true"}});
  model3.precompute(model_params);
  EXPECT_NEAR(0.037037403706123712, model3.energy(site.position(), site, *config, model_params), 1e-15);

  pos.set_vector({0, 0, -3.});
  config->displace_particle(config->selection_of_all(), pos);
  EXPECT_LT(1e100, model.energy(site.position(), site, *config, model_params));
}

}  // namespace feasst

#include <gtest/gtest.h>
#include "confinement/include/model_hard_shape.h"
#include "core/include/configuration.h"
#include "core/include/select_list.h"

namespace feasst {

TEST(ModelHardShape, half_space) {
  auto half_space = HalfSpace()
    .set_dimension(2)
    .set_intersection(1)
    .set_direction(1);
  ModelHardShape model;
  model.set_shape(std::make_shared<HalfSpace>(half_space));
  Configuration config;
  config.set_domain(Domain().set_cubic(8));
  config.add_particle_type("../forcefield/data.atom");
  config.add_particle_of_type(0);
  const ModelParams model_params = config.unique_types().model_params();
  EXPECT_LT(1e100, model.energy(config.particle(0).site(0), &config, model_params));
  Position pos;
  pos.set_vector({-24.23, 35.45, 1.5000001});
  config.displace_particles(config.selection_of_all(), pos);
  EXPECT_NEAR(0, model.energy(config.particle(0).site(0), &config, model_params), 1e-15);
  pos.set_vector({0, 0, -.00001});
  config.displace_particle(config.selection_of_all(), pos);
  EXPECT_LT(1e100, model.energy(config.particle(0).site(0), &config, model_params));
}

}  // namespace feasst

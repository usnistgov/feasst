#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "system/include/potential.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "morph/include/perturb_particle_type.h"
#include "charge/include/utils.h"
#include "charge/test/charge_utils.h"

namespace feasst {

TEST(PerturbParticleType, serialize) {
  System sys;
  {
    auto config = MakeConfiguration({{"cubic_side_length", "8"},
      {"particle_type0", "../particle/lj.txt"}});
    config->add_particle_type("../particle/lj.txt", "sig0.25");
    config->add_particle_of_type(0);
    config->add_particle_of_type(0);
    config->set_model_param("sigma", 1, 0.25);
    config->set_model_param("cutoff", 1, 1.0);
    config->update_positions({{0, 0, 0}, {1, 1, 1}});
    config->particle_type_to_group_create(0);
    config->particle_type_to_group_create(1);
    EXPECT_EQ(config->model_params().select("sigma").mixed_value(0, 1), 1.25/2.);
    EXPECT_EQ(config->model_params().select("cutoff").mixed_value(0, 1), 2.);
    sys.add(config);
  }
  sys.add(MakePotential(MakeLennardJones()));
  const Configuration& config = sys.configuration();
  auto morph = MakePerturbParticleType({{"type", "1"}});
  PerturbParticleType morph2 = test_serialize(*morph);
  EXPECT_EQ(config.particle(0).type(), 0);
  EXPECT_EQ(config.particle(0).site(0).type(), 0);
  EXPECT_EQ(config.num_particles_of_type(0), 2);
  EXPECT_EQ(config.num_particles_of_type(1), 0);
  EXPECT_EQ(config.group_select(0).num_particles(), 2);
  EXPECT_EQ(config.group_select(1).num_particles(), 2);
  EXPECT_EQ(config.group_select(2).num_particles(), 0);
  auto random = MakeRandomMT19937();
  auto select = MakeTrialSelectParticle({{"particle_type", "0"}});
  select->precompute(&sys);
  select->sel(&sys, random.get());
  EXPECT_NEAR(sys.energy(), -0.142661179698217000, NEAR_ZERO);
  morph2.perturb(&sys, select.get(), random.get());
  sys.configuration().check();
  const Particle& part = config.select_particle(select->mobile().particle_index(0));
  EXPECT_EQ(part.type(), 1);
  EXPECT_EQ(part.site(0).type(), 1);
  EXPECT_EQ(config.num_particles_of_type(0), 1);
  EXPECT_EQ(config.num_particles_of_type(1), 1);
  EXPECT_EQ(config.group_select(0).num_particles(), 2);
  EXPECT_EQ(config.group_select(1).num_particles(), 1);
  EXPECT_EQ(config.group_select(2).num_particles(), 1);
  EXPECT_NEAR(sys.energy(), -0.0088108241166350975, NEAR_ZERO);
  morph2.revert(&sys);
  EXPECT_EQ(config.particle(0).type(), 0);
  EXPECT_EQ(config.particle(0).site(0).type(), 0);
  EXPECT_EQ(config.num_particles_of_type(0), 2);
  EXPECT_EQ(config.num_particles_of_type(1), 0);
  EXPECT_EQ(config.group_select(0).num_particles(), 2);
  EXPECT_EQ(config.group_select(1).num_particles(), 2);
  EXPECT_EQ(config.group_select(2).num_particles(), 0);
  EXPECT_NEAR(sys.energy(), -0.142661179698217000, NEAR_ZERO);
}

// likely unnecessary test, could delete this
TEST(PerturbParticleType, ewald) {
  System sys = rpm({{"alpha", str(5.6/12)}, {"kmax_squared", "38"}});
  sys.get_configuration()->add_particle_of_type(0);
  sys.get_configuration()->add_particle_of_type(1);
  sys.get_configuration()->update_positions({{0, 0, 0}, {1, 1, 1}});
  const double en_init = sys.energy();
  EXPECT_NEAR(en_init, -0.58032644318603566, NEAR_ZERO);
  EXPECT_EQ(sys.configuration().num_particles_of_type(0), 1);
  EXPECT_EQ(sys.configuration().num_particles_of_type(1), 1);

  auto select = MakeTrialSelectParticle({{"particle_type", "0"}});
  select->set_trial_state(0);
  select->precompute(&sys);
  auto random = MakeRandomMT19937();
  select->sel(&sys, random.get());
  auto morph = MakePerturbParticleType({{"type", "1"}});
  const double en_old = sys.perturbed_energy(select->mobile());
  morph->perturb(&sys, select.get(), random.get());
  select->set_trial_state(1);
  const double en_new = sys.perturbed_energy(select->mobile());
  const double en_final = sys.energy();
  EXPECT_NEAR(en_init + en_new - en_old, en_final, NEAR_ZERO);

  EXPECT_NEAR(en_final, 0.29630027798728265, NEAR_ZERO);
  EXPECT_EQ(sys.configuration().num_particles_of_type(0), 0);
  EXPECT_EQ(sys.configuration().num_particles_of_type(1), 2);
}

}  // namespace feasst

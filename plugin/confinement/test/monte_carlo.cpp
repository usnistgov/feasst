#include "utils/test/utils.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "confinement/include/sphere.h"
#include "confinement/include/slab.h"
#include "confinement/include/shape_union.h"
#include "confinement/include/cylinder.h"
#include "confinement/include/model_hard_shape.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "steppers/include/movie.h"
#include "steppers/include/log.h"

namespace feasst {

TEST(MonteCarlo, ShapeUnion) {
  MonteCarlo mc;
  mc.add(Configuration(MakeDomain({{"cubic_box_length", "8"}}),
    {{"particle_type", "../forcefield/data.lj"}}));
  mc.add(Potential(MakeLennardJones()));
  mc.add(Potential(MakeModelHardShape(MakeShapeUnion(
    MakeSphere(
      {{"radius", "2"}},
      Position({{"x", "0"}, {"y", "0"}, {"z", "0"}})),
    MakeSlab({{"dimension", "2"}, {"bound0", "-1"}, {"bound1", "1"}})
  ))));
  mc.add(MakeMetropolis({
    {"beta", "1.2"},
    {"chemical_potential", "1."},
  }));
  mc.add(MakeTrialTranslate({
    {"weight", "1."},
    {"tunable_param", "2."},
  }));
  SeekNumParticles(10).with_trial_add().run(&mc);
  const int steps_per = 1e0;
  mc.add(MakeLog({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/confine_log.txt"}
  }));
  mc.add(MakeMovie({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/confine_movie.xyz"}
  }));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(1e3);
}

std::shared_ptr<Shape> porous_network() {
  auto network = MakeShapeUnion(
    MakeSphere(
      {{"radius", "5"}},
      Position({{"x", "0"}, {"y", "0"}, {"z", "0"}})),
    MakeCylinder({{"radius", "2"}},
      Position({{"x", "0"}, {"y", "0"}, {"z", "0"}}),
      Position({{"x", "0"}, {"y", "0"}, {"z", "1"}})));
  network = MakeShapeUnion(network,
    MakeCylinder({{"radius", "2"}},
      Position({{"x", "0"}, {"y", "0"}, {"z", "0"}}),
      Position({{"x", "0"}, {"y", "1"}, {"z", "0"}})));
  network = MakeShapeUnion(network,
    MakeCylinder({{"radius", "2"}},
      Position({{"x", "0"}, {"y", "0"}, {"z", "0"}}),
      Position({{"x", "1"}, {"y", "0"}, {"z", "0"}})));
  return network;
}

TEST(MonteCarlo, ShapeUnion_LONG) {
  MonteCarlo mc;
  mc.add(Configuration(MakeDomain({{"cubic_box_length", "20"}}),
    {{"particle_type", "../forcefield/data.lj"}}));
  mc.add(Potential(MakeLennardJones()));
  mc.add(Potential(MakeModelHardShape(porous_network())));
  mc.add(MakeMetropolis({
    {"beta", "1.5"},
    {"chemical_potential", "1."},
  }));
  mc.add(MakeTrialTranslate({
    {"weight", "1."},
    {"tunable_param", "2."},
  }));
  SeekNumParticles(500).with_trial_add().run(&mc);
  const int steps_per = 1e5;
  mc.add(MakeLog({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/confine_log.txt"}
  }));
  mc.add(MakeMovie({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/confine_movie.xyz"},
    {"clear_file", "true"},
  }));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(1e6);
}

}  // namespace feasst

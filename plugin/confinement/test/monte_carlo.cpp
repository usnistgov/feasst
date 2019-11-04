#include "utils/test/utils.h"
#include "confinement/include/sphere.h"
#include "confinement/include/slab.h"
#include "confinement/include/model_hard_shape.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trial_translate.h"
#include "steppers/include/movie.h"
#include "steppers/include/log.h"
#include "system/include/model_lj.h"
#include "monte_carlo/include/metropolis.h"

namespace feasst {

TEST(MonteCarlo, ShapeUnion) {
  MonteCarlo mc;
  mc.add(Configuration({
    {"cubic_box_length", "8"},
    {"particle_type", "../forcefield/data.lj"},
  }));
  mc.add(Potential(MakeModelLJ()));
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
  mc.seek_num_particles(10);
  const int steps_per = 1e0;
  mc.add(MakeLog({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/confine_log.txt"}
  }));
  mc.add(MakeMovie({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/confine_movie.xyz"}
  }));
  mc.attempt(1e3);

  // serialize
  test_serialize(mc);
}

}  // namespace feasst

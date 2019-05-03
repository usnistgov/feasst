
#include <gtest/gtest.h>
#include "confinement/include/sphere.h"
#include "confinement/include/slab.h"
#include "confinement/include/model_hard_shape.h"
#include "core/include/monte_carlo.h"
#include "core/include/trial_translate.h"
#include "core/include/movie.h"
#include "core/include/log.h"
#include "core/include/model_lj.h"
#include "core/include/criteria_metropolis.h"

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
  mc.add(MakeCriteriaMetropolis({
    {"beta", "1.2"},
    {"chemical_potential", "1."},
  }));
  mc.add(MakeTrialTranslate({
    {"weight", "1."},
    {"max_move", "2."},
  }));
  mc.seek_num_particles(50);
  const int steps_per = 1e3;
  mc.add(MakeLog({{"steps_per", str(steps_per)}}));
  mc.add(MakeMovie(
   {{"steps_per", str(steps_per)},
    {"file_name", "movie.xyz"}}));
  mc.attempt(1e3);
}

}  // namespace feasst

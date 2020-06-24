#include "utils/test/utils.h"
#include "utils/include/checkpoint.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/configuration.h"
#include "configuration/include/select.h"
#include "configuration/include/domain.h"
#include "shape/include/slab_sine.h"
#include "confinement/include/model_table_cartesian.h"

namespace feasst {

TEST(ModelTableCart3FoldSym, SineSlabTable_LONG) {
  auto hamaker = MakeModelTableCart3FoldSym(MakeTable3D({
    {"num0", "11"},
    {"num1", "11"},
    {"num2", "11"},
    {"default_value", "0."}}));
  #ifdef _OPENMP
  // HWH something wrong with bisection and OMP
  hamaker->compute_table(
  #elif // _OPENMP
  hamaker->compute_table(
  #endif // _OPENMP
    MakeSlabSine(MakeFormulaSineWave({{"amplitude", "2"}, {"width", "8"}}),
      { {"dimension", "0"}, {"wave_dimension", "1"}, {"average_bound0", "-5"},
        {"average_bound1", "5"}}).get(),
    MakeDomain({{"cubic_box_length", "20"}}).get(),
    MakeRandomMT19937({{"seed", "123"}}).get(),
    { {"alpha0", "6"},
      {"epsilon0", "-1"},
      {"alpha1", "12"},
      {"epsilon1", "1"},
      {"max_radius", "10"},
      {"num_shells", "100"},
      {"points_per_shell", "100"}});
  MakeCheckpoint({{"file_name", "tmp/sine_slab_table"}})->write(hamaker->table());
}

}  // namespace feasst

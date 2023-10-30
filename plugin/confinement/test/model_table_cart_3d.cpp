#include "utils/test/utils.h"
#include "utils/include/checkpoint.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/configuration.h"
#include "configuration/include/select.h"
#include "configuration/include/domain.h"
#include "shape/include/slab_sine.h"
#include "shape/include/formula_sine_wave.h"
#include "confinement/include/model_table_cartesian.h"

namespace feasst {

TEST(ModelTableCart3DIntegr, SineSlabTable_LONG) {
  auto domain = MakeDomain({{"cubic_side_length", "20"}});
  auto hamaker = MakeModelTableCart3DIntegr(MakeTable3D({
    {"num0", "3"},
    {"num1", "3"},
    {"num2", "3"},
    {"default_value", "0."}}));
  #ifdef _OPENMP
  // HWH something wrong with bisection and OMP
  hamaker->compute_table(
  #else // _OPENMP
  hamaker->compute_table(
  #endif // _OPENMP
    MakeSlabSine({ {"amplitude", "2"}, {"width", "8"},
      {"dimension", "0"}, {"wave_dimension", "1"}, {"average_bound0", "-5"},
      {"average_bound1", "5"}}).get(),
    *domain,
    MakeRandomMT19937({{"seed", "123"}}).get(),
    { {"alpha0", "6"},
      {"epsilon0", "-1"},
      {"alpha1", "12"},
      {"epsilon1", "1"},
      {"max_radius", "10"},
      {"num_shells", "10"},
      {"points_per_shell", "10"}});
  MakeCheckpoint({{"checkpoint_file", "tmp/sine_slab_table"}})->write(hamaker->table());
  ModelTableCart3DIntegr hamaker2 = test_serialize(*hamaker);
  hamaker2.write("tmp/table3d_out.txt");
  auto hamaker3 = MakeModelTableCart3DIntegr({{"table_file", "tmp/table3d_out.txt"}});
}

TEST(ModelTableCart3DIntegr, SineSlabTable_TXT_LONG) {
//  MakeDomain({{"cubic_side_length", "20"}}).get(),
  auto hamaker = MakeModelTableCart3DIntegr({
    {"table_file", "tmp/table3d_out.txt"},
    {"cubic_side_length", "20"},
    {"num0", "11"},
    {"num1", "11"},
    {"num2", "11"},
    {"default_value", "0."},
  #ifdef _OPENMP
    {"use_omp", "true"},
  #else // _OPENMP
    {"use_omp", "false"},
  #endif // _OPENMP
    {"shape_file", "../plugin/confinement/test/data/shape.txt"},
    {"random", "RandomMT19937"},
    {"alpha0", "6"},
    {"epsilon0", "-1"},
    {"alpha1", "12"},
    {"epsilon1", "1"},
    {"max_radius", "10"},
    {"num_shells", "10"},
    {"points_per_shell", "10"}});
  auto config = MakeConfiguration({{"cubic_side_length", "20"}});
  hamaker->precompute(config->model_params());
  MakeCheckpoint({{"checkpoint_file", "tmp/sine_slab_table"}})->write(hamaker->table());
  ModelTableCart3DIntegr hamaker2 = test_serialize(*hamaker);
  hamaker2.write("tmp/table3d_out.txt");
  auto hamaker3 = MakeModelTableCart3DIntegr({{"table_file", "tmp/table3d_out.txt"}});
}

TEST(ModelTableCard1DHard, compute_table) {
  auto domain = MakeDomain({{"cubic_side_length", "10"}});
  auto shape = MakeSlabSine(
    { {"dimension", "1"}, {"wave_dimension", "0"}, {"average_bound0", "-5"},
      {"average_bound1", "5"}, {"amplitude", "2"}, {"width", "10"}});
  auto table = MakeTable1D({{"num", "11"}});
  auto model = MakeModelTableCart1DHard(table);
  auto random = MakeRandomMT19937();
  model->compute_table(shape.get(), domain.get(), random.get());
}

}  // namespace feasst

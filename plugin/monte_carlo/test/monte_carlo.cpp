#include <memory>
#include "utils/test/utils.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "monte_carlo/include/criteria_metropolis.h"
#include "utils/include/utils_io.h"
#include "math/include/accumulator.h"
#include "system/test/system_test.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"

namespace feasst {

TEST(MonteCarlo, serialize) {
  MonteCarlo mc = mc_lj(), mc2 = test_serialize(mc);
}

TEST(MonteCarlo, NVT_benchmark) {
  seed_random_by_date();
  MonteCarlo mc = mc_lj();
  mc.seek_num_particles(50);
  // mc.attempt(1e6);  // ~4 seconds
  mc.attempt(1e4);
}

TEST(MonteCarlo, NVT_SRSW) {
  seed_random_by_date();
  MonteCarlo mc = mc_lj();
  const int nMol = 500;
  const double rho = 1e-3, length = pow(static_cast<double>(nMol)/rho, 1./3.);
  mc.get_system()->get_configuration()->set_side_length(
    Position().set_vector({length, length, length}));
  mc.seek_num_particles(nMol);
  Accumulator pe;
  for (int trial = 0; trial < 1e3; ++trial) {
    mc.attempt(1);  // ~4 seconds
    pe.accumulate(mc.criteria()->current_energy());
  }
  // HWH temperature not set
  INFO("pe " << pe.average());
}

TEST(MonteCarlo, GCMC) {
  MonteCarlo mc = mc_lj();
  mc.add(MakeTrialTransfer());
  mc.attempt(1e4);
}

}  // namespace feasst

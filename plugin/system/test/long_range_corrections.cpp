#include <cmath>
#include "utils/test/utils.h"
#include "configuration/include/domain.h"
#include "configuration/test/config_utils.h"
#include "system/include/long_range_corrections.h"
#include "system/include/model_empty.h"

namespace feasst {

TEST(LongRangeCorrections, LRC) {
  Configuration config = two_particle_configuration();
  ModelEmpty empty;
  LongRangeCorrections lrc;
  lrc.precompute(&config);
  empty.compute(&config, &lrc);
  const double pe_lrc = (8./3.)*PI*std::pow(config.num_particles(), 2)/config.domain().volume()
    *((1./3.)*std::pow(3, -9) - std::pow(3, -3));
  EXPECT_NEAR(pe_lrc, lrc.energy(), NEAR_ZERO);

  auto visit = test_serialize<LongRangeCorrections, VisitModel>(lrc);
  empty.compute(&config, visit.get());
  EXPECT_NEAR(visit->energy(), lrc.energy(), NEAR_ZERO);
}

TEST(LongRangeCorrections, many_types) {
  auto config = MakeConfiguration({{"cubic_side_length", "8.5"},
    {"particle_type0", "../particle/atom.txt"},
    {"particle_type1", "../particle/lj.txt"},
    {"add_particles_of_type0", "4"},
    {"add_particles_of_type1", "3"}});
  ModelEmpty empty;
  LongRangeCorrections lrc;
  lrc.precompute(config.get());
  empty.compute(config.get(), &lrc);
  const double pe_lrc = (8./3.)*PI*std::pow(config->num_particles(), 2)/config->domain().volume()
    *((1./3.)*std::pow(3, -9) - std::pow(3, -3));
  EXPECT_NEAR(pe_lrc, lrc.energy(), NEAR_ZERO);

  auto visit = test_serialize<LongRangeCorrections, VisitModel>(lrc);
  empty.compute(config.get(), visit.get());
  EXPECT_NEAR(visit->energy(), lrc.energy(), NEAR_ZERO);
}

}  // namespace feasst

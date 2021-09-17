#include <cmath>
#include "utils/test/utils.h"
#include "configuration/include/domain.h"
#include "configuration/include/utils.h"
#include "system/include/long_range_corrections.h"
#include "system/include/model_empty.h"

namespace feasst {

TEST(LongRangeCorrections, LRC) {
  Configuration config = two_particle_configuration();
  ModelEmpty empty;
  LongRangeCorrections lrc;
  empty.compute(&config, &lrc);
  const double pe_lrc = (8./3.)*PI*std::pow(config.num_particles(), 2)/config.domain().volume()
    *((1./3.)*std::pow(3, -9) - std::pow(3, -3));
  EXPECT_NEAR(pe_lrc, lrc.energy(), NEAR_ZERO);

  auto visit = test_serialize<LongRangeCorrections, VisitModel>(lrc);
  empty.compute(&config, visit.get());
  EXPECT_NEAR(visit->energy(), lrc.energy(), NEAR_ZERO);
}

TEST(LongRangeCorrections, many_types) {
  Configuration config(MakeDomain({{"cubic_box_length", "8.5"}}),
    {{"particle_type0", "../forcefield/atom.fstprt"},
     {"particle_type1", "../forcefield/lj.fstprt"}});
  for (int i = 0; i < 4; ++i) config.add_particle_of_type(0);
  for (int i = 0; i < 3; ++i) config.add_particle_of_type(1);
  ModelEmpty empty;
  LongRangeCorrections lrc;
  empty.compute(&config, &lrc);
  const double pe_lrc = (8./3.)*PI*std::pow(config.num_particles(), 2)/config.domain().volume()
    *((1./3.)*std::pow(3, -9) - std::pow(3, -3));
  EXPECT_NEAR(pe_lrc, lrc.energy(), NEAR_ZERO);

  auto visit = test_serialize<LongRangeCorrections, VisitModel>(lrc);
  empty.compute(&config, visit.get());
  EXPECT_NEAR(visit->energy(), lrc.energy(), NEAR_ZERO);
}

}  // namespace feasst

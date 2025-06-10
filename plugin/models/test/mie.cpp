#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "system/include/long_range_corrections.h"
#include "system/include/model_empty.h"
#include "models/include/mie.h"

namespace feasst {

TEST(Mie, analytical) {
  auto config = MakeConfiguration({{"particle_type0", "../particle/mie.txt"},
                                   {"add_particles_of_type0", "2"},
                                   {"cubic_side_length", "8"}});
  auto model1 = MakeMie();
  model1->precompute(*config);
  std::shared_ptr<Model> model2 = test_serialize<Mie, Model>(*model1, "Mie 2094 1 0 2 -1 2905 3 4 mie_prefactor 4795 1 0 1 1 4.9207071226910948 1 1 1 7964 ");
  DEBUG(model2->energy(1.5*1.5, 0, 0, config->model_params()));
  EXPECT_NEAR(-0.17514250679168769, model2->energy(1.5*1.5, 0, 0, config->model_params()), NEAR_ZERO);

  // test LRCs
  auto lrc = MakeLongRangeCorrections();
  lrc->precompute(config.get());
  auto lrc2 = test_serialize<LongRangeCorrections, VisitModel>(*lrc);
  ModelEmpty empty;
  empty.compute(config.get(), lrc2.get());
  EXPECT_NEAR(-0.000198678220854467, lrc2->energy(), NEAR_ZERO);
}

}  // namespace feasst

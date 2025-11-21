#include "utils/test/utils.h"
#include "gibbs/include/trial_gibbs_morph.h"

namespace feasst {

TEST(TrialGibbsMorphOneWay, serialize) {
  auto add = std::make_shared<TrialGibbsMorphOneWay>(argtype({{"to_config", "0"}, {"particle_type", "0"}, {"particle_type_morph", "1"}}));
  std::shared_ptr<Trial> trial2 = test_serialize<TrialGibbsMorphOneWay, Trial>(*add);
}

}  // namespace feasst

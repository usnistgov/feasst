#include "utils/test/utils.h"
#include "gibbs/include/trial_gibbs_particle_transfer.h"

namespace feasst {

TEST(TrialGibbsParticleTransferOneWay, serialize) {
  auto add = MakeTrialGibbsParticleTransferOneWay({{"to_config", "0"}});
  std::shared_ptr<Trial> trial2 = test_serialize<TrialGibbsParticleTransferOneWay, Trial>(*add);
}

}  // namespace feasst

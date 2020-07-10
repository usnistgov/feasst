#include "utils/test/utils.h"
#include "cluster/include/trial_avb4.h"

namespace feasst {

TEST(TrialAVB4, serialize) {
  TrialAVB4 trial(MakeNeighborCriteria({{"maximum_distance", "3"},
                                        {"minimum_distance", "1"}}));
  TrialAVB4 trial2 = test_serialize(trial);
}

}  // namespace feasst

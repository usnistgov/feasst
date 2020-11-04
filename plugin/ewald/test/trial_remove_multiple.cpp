#include "utils/test/utils.h"
#include "ewald/include/trial_remove_multiple.h"

namespace feasst {

TEST(TrialRemoveMultiple, serialize) {
  auto remove = MakeTrialRemoveMultiple({{"particle_type0", "2"},
                               {"particle_type1", "3"}});
  EXPECT_EQ(remove->stage(0).trial_select().particle_type(), 2);
  EXPECT_EQ(remove->stage(1).trial_select().particle_type(), 3);
//  std::stringstream ss;
//  remove.serialize(ss);
//  INFO(ss.str());
//  TrialRemoveMultiple remove2(ss);
//  remove2.serialize(ss);
//  INFO(ss.str());
  Trial remove2 = test_serialize(*remove);
}

}  // namespace feasst

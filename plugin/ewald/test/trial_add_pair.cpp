#include "utils/test/utils.h"
#include "ewald/include/trial_add_pair.h"

namespace feasst {

TEST(TrialAddPair, serialize) {
  auto add = MakeTrialAddPair({{"particle_type0", "2"},
                               {"particle_type1", "3"}});
  EXPECT_EQ(add->stage(0)->trial_select()->particle_type(), 2);
  EXPECT_EQ(add->stage(1)->trial_select()->particle_type(), 3);
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialAddPair add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialAddPair add2 = test_serialize(*add);
}

}  // namespace feasst

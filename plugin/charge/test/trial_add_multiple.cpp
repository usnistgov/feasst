#include "utils/test/utils.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_select.h"
#include "charge/include/trial_add_multiple.h"

namespace feasst {

TEST(TrialAddMultiple, serialize) {
  auto add = MakeTrialAddMultiple({{"particle_type0", "2"},
                                   {"particle_type1", "3"},
                                   {"reference_index", "-1"}});
  EXPECT_EQ(add->stage(0).trial_select().particle_type(), 2);
  EXPECT_EQ(add->stage(1).trial_select().particle_type(), 3);
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialAddMultiple add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  Trial add2 = test_serialize(*add);
}

}  // namespace feasst

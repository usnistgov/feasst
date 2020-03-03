#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "chain/include/trial_select_reptate.h"
#include "chain/test/system_chain.h"

namespace feasst {

TEST(TrialSelectReptate, serialize) {
  auto add = MakeTrialSelectReptate({{"max_length", "1"}});
  auto sys = chain_system();
  add->precompute(&sys);
  auto random = MakeRandomMT19937();
  add->select(Select(), &sys, random.get());
  DEBUG(add->mobile().str());
  std::stringstream ss;
  add->serialize(ss);
  DEBUG(ss.str());
//  TrialSelectReptate add2(ss);
//  add2.serialize(ss);
//  DEBUG(ss.str());
  TrialSelectReptate add2 = test_serialize(*add);
  auto sel = test_serialize<TrialSelectReptate, TrialSelect>(*add);
  //std::shared_ptr<TrialSelect>
}

}  // namespace feasst

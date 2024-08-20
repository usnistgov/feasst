#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "system/include/potential.h"
#include "system/include/system.h"
#include "chain/include/select_reptate.h"
#include "chain/test/system_chain.h"

namespace feasst {

TEST(SelectReptate, serialize) {
  auto add = MakeSelectReptate({{"max_length", "1"}, {"particle_type", "0"}});
  auto sys = chain_system();
  add->precompute(&sys);
  auto random = MakeRandomMT19937();
  add->sel(&sys, random.get());
  DEBUG(add->mobile().str());
  std::stringstream ss;
  add->serialize(ss);
  DEBUG(ss.str());
//  SelectReptate add2(ss);
//  add2.serialize(ss);
//  DEBUG(ss.str());
  SelectReptate add2 = test_serialize(*add);
  auto sel = test_serialize<SelectReptate, TrialSelect>(*add);
  //std::shared_ptr<Select>
}

}  // namespace feasst

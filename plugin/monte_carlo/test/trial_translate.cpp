#include "utils/test/utils.h"
#include "monte_carlo/include/trial_translate.h"

namespace feasst {

TEST(TrialTranslate, serialize) {
  TrialTranslate add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialTranslate add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialTranslate add2 = test_serialize(add);
}

}  // namespace feasst
